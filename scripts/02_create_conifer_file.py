import numpy as np
import math
import operator
import random
import time
import sys
import itertools
import re
import os
from tables import *
import argparse
import logging
import pandas 

def process_sample_list(args):
    if os.path.isdir(args.samples):
        # got a directory of samples, hopefully RPKM files
        print "Not yet supported!"
        sys.exit(1)
    else:
        try:
            with open(args.samples) as sample_file:
                # read in CSV file. Format should be SampleID, RPKM_PATH, all fields after are additional fields added.
                sample_list = pandas.read_csv(sample_file)
        except:
            log.exception("Could not read sample file!")
            sys.exit(1)
    return sample_list

def convert_chrom(x):
    t = str(x).lstrip("chr")
    if t == "X":
        return 23
    elif t == "Y":
        return 24
    else:
        return int(t)

def load_probes(probes, minsize=None):
    """
    Load and format a probe file, optionally expanding small probes to <minsize>.
    """
    try:
        probes = pd.read_csv(probes, sep="\t", names=["chromosome", "start", "stop", "name", "isSegDup","isPPG","strand"])
    except:
        log.exception("Could not read probes file. "
                        "Check file exists, is tab-delimited and has appropriate header?"
                        "Required fields: chromosome,start,stop,name,isSegDup,isPPG,strand")

    probes["chromosome"] = map(convert_chrom, probes.chromosome)
    if probes["chromosome"].dtype != np.int64:
        log.exception("Could not convert all probes to standard chromosome notation!")

    probes["probe_size"] = probes["stop"] - probes["start"]
    probes["probeID"] = probes.index
    if minsize:
        mask = probes.probe_size < minsize
        d = minsize - probes[mask]["probe_size"]
        left_pad = np.floor(d/2)
        right_pad = np.ceil(d/2)
        probes["start"][mask] = probes["start"][mask] - left_pad
        probes["stop"][mask] = probes["stop"][mask] + right_pad

    return probes

def zrpkm(rpkm,median,sd):
    return (rpkm - median) / sd

class rpkm_value(IsDescription):
    probeID = UInt32Col(pos=0)
    rpkm = FloatCol(pos=1)

class probe(IsDescription):
    probeID =   UInt32Col(pos=0)
    start =     UInt32Col(pos=1) # start of probe
    stop =      UInt32Col(pos=2) # stop of probe
    name =      StringCol(20,pos=3)   # 20-character String
    isSegDup =  BoolCol(pos=4)  # is probe in segdup
    isPPG  =    BoolCol(pos=5)  # is gene a processed pseudogene
    strand =    StringCol(1,pos=6)  # strand of gene

class sample(IsDescription):
    sampleID =  StringCol(20,pos=0) # 20-char string (sampleID)
    cohort =    StringCol(20,pos=1)
    sex =       StringCol(1,pos=2)
    ethnicity = StringCol(50,pos=3)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--outfile", "-O", type=str, required=True)
    parser.add_argument("--components_removed", "-C", type=int, required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--probes", required=True)
    parser.add_argument("--chromosomes",type=int, nargs="*", required=False, default=range(1,25))
    parser.add_argument("--min-probe-size", default=10, type=int)
    parser.add_argument("--QC_report","--qc", type=str, required=False, default=None)
    parser.add_argument("--log", type=str, required=False, default=None)
    parser.add_argument("--loglevel", type=str, required=False, default="INFO")
    args = parser.parse_args()
    


    numeric_log_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    try:
        logging.basicConfig(filename=args.log,level=numeric_log_level,format='[%(levelname)s] [%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    except IOError:
        print "Could not open log file or start logging. Log output will be sent to stderr/stdout"
        logging.basicConfig(level=numeric_log_level,format='[%(levelname)s] [%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    log = logging.getLogger("CoNIFER")
    log.debug(args)

    try:
        h5file_out = openFile(args.outfile, mode='w')
    except:
        log.exception("Could not open CoNIFER outfile (%s) for writing" % args.out)
        sys.exit(1)


    if args.QC_report is not None:
        qc_report = True
        try:
            QC_report_file=open(args.QC_report,"w")
        except:
            log.exception("Could not open QC report file (%s) for writing" % args.QC_report)
            sys.exit(1)
    else:
        qc_report = False

    log.info("Processing sample list...")
    sample_list = process_sample_list(args)
    num_samples = len(sample_list)
    log.info("Have %d samples", len(sample_list))

    log.info("Processing probe list...")
    probe_list = load_probes(args.probes, minsize=args.min_probe_size)
    num_probes = len(probe_list)
    log.info("Have %d probes", num_probes)

    probe_group = h5file_out.createGroup("/","probes","probes")

    for chromosome in args.chromosomes:
        log.info("Processing chromosome %d" % chromosome)

        chr_probe_mask = probe_list.chromosome == chromosome
        chr_num_probes = np.sum(chr_probe_mask)
        log.debug("Chromosome %d has %d probes", chromosome, chr_num_probes)

        chr_group_name = "chr%d" % chromosome
        chr_group = h5file_out.createGroup("/",chr_group_name,chr_group_name)

        start_probeID = min(probe_list[chr_probe_mask]['probeID'].values)
        stop_probeID = max(probe_list[chr_probe_mask]['probeID'].values)
        log.debug("Chromosome %d  probeID range: %d - %d", chromosome, start_probeID, stop_probeID)
        
        rpkm = np.zeros([chr_num_probes,num_samples], dtype=np.float)

        log.info("Reading RPKM files")
        cnt = 0
        out_samples = []
        for ix,s in sample_list.iterrows():
            rpkm_infile_fn = s["file"]
            try:
                h5file_in = openFile(rpkm_infile_fn, mode='r')
                #log.debug("Loading %s (%d/%d)" % (rpkm_infile_fn, ix+1, num_samples))
            except IOError: 
                #print "cannot find %s.h5" % s
                log.warn("Cannot open RPKM data file: %s. Removing from analysis.", rpkm_infile_fn)
                continue
            #print s
            out_samples.append(str(s["sampleID"]))
            d = h5file_in.root.rpkm.rpkm.read(field="rpkm",start=start_probeID-1,stop=stop_probeID)
            if np.sum(np.isnan(d)) > 0:
                log.warn("Sample %s contains %d null values! Removing from analysis" % (s["sampleID"], np.sum(np.isnan(d))))
                continue

            rpkm[:,cnt] = d
            cnt +=1
            h5file_in.close()
        
        log.info("Successfully read %d of %d input RPKM files", cnt, num_samples)
        #shrink rpkm array to only size of input (for failed files)
        try:
            rpkm = rpkm[:,0:len(out_samples)]
        except MemoryError:
            log.exception("Ran out of memory in allocating RPKM data matrix!")
            sys.exit(1)
        
        # get the SD and median
        median = np.median(rpkm,1)
        sd = np.std(rpkm,1)
        
        # filter by median rpkm >= 1
        median_probe_mask = median >= 1
        rpkm = rpkm[median_probe_mask, :]
        num_masked_probes = np.sum(median_probe_mask)
        log.info("Masking %d probes of %d", chr_num_probes-num_masked_probes, chr_num_probes)

        dt = np.dtype([('probeID',np.uint32),('start',np.uint32),('stop',np.uint32), ('name', np.str_, 20),('isSegDup', np.bool),('isPPG', np.bool),('strand', np.str_,1)])
        out_probes = np.empty(num_masked_probes,dtype=dt)

        for field in ["probeID", "start", "stop", "name", "isSegDup", "isPPG", "strand"]:
            out_probes[field] = probe_list[chr_probe_mask][median_probe_mask][field].values
        
        probe_table = h5file_out.createTable(probe_group,"probes_chr%d" % chromosome,probe,"chr%d" % chromosome)
        probe_table.append(out_probes)
        
        
        # Z-score using zrpkm function above
        log.info("Z-transforming data")
        rpkm = np.apply_along_axis(zrpkm, 0, rpkm, median[median_probe_mask], sd[median_probe_mask])
        
        # svd transform
        log.info("SVD transform of data")
        try:
            U, S, Vt = np.linalg.svd(rpkm,full_matrices=False)
        except MemoryError:
            log.exception("Ran out of memory in SVD transformation!")
            sys.exit(1)
        except np.linalg.linalg.LinAlgError:
            print rpkm.shape
            np.save("/net/eichler/vol8/home/nkrumm/zrpkm.tmp.svd", rpkm)

        new_S = np.diag(np.hstack([np.zeros([args.components_removed]),S[args.components_removed:]]))
        
        if qc_report:
            QC_report_file.write('chr' + str(chromosome) + '\t' + '\t'.join([str(_i) for _i in S]) + "\n")

        log.info("Reconstructing original data matrix")
        # reconstruct data matrix
        try:
            rpkm = np.dot(U, np.dot(new_S, Vt))
        except MemoryError:
            log.exception("Ran out of memory in reconstructing data matrix!")
            sys.exit(1)

        log.info("Done with CoNIFER transformation")
        
        # save to HDF5 file
        log.info("Saving values to file")
        probeIDs  = probe_list[chr_probe_mask][median_probe_mask]["probeID"].values
        for i,sampleID in enumerate(out_samples):
            out_data = np.empty(num_masked_probes,dtype='u4,f8')
            out_data['f0'] = probeIDs
            out_data['f1'] = rpkm[:,i]
            sample_tbl = h5file_out.createTable(chr_group,"sample_" + str(sampleID),rpkm_value,"%s" % sampleID)
            sample_tbl.append(out_data)

        log.info("Finished with chromosome %d", chromosome)


log.info("Adding sample table information to file")
sample_group = h5file_out.createGroup("/","samples","samples")
sample_table = h5file_out.createTable(sample_group,"samples",sample,"samples")

dt = np.dtype([('sampleID',np.str_,20),('cohort',np.str_,20),('sex',np.str_,1), ('ethnicity', np.str_, 50)])
# create mask of sample actually in HDF5 file:
mask = np.array([str(s) in out_samples for s in sample_list.sampleID.values])

out_samples = np.empty(np.sum(mask),dtype=dt)

for field in ["sampleID", "cohort", "sex", "ethnicity"]:
    out_samples[field] = sample_list[mask][field]

sample_table.append(out_samples)


log.info("Closing HDF5 file and other open files")
h5file_out.close()
if qc_report:
    QC_report_file.close()
log.info("Finished.")
