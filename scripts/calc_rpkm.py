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
import pandas as pd
import argparse

class exon(IsDescription):
    exonID = UInt32Col(pos=0)
    read_cnt = UInt32Col(pos=1)
    rpkm = FloatCol(pos=2)

def convert_chrom(x):
    t = x.lstrip("chr")
    if t == "X":
        return 23
    elif t == "Y":
        return 24
    else:
        return int(t)

def load_probes(probe_file, minsize=None):
    """
    Load and format a probe file, optionally expanding small probes to <minsize>.
    """
    try:
        probes = pd.read_csv(probe_file, sep="\t", names=["chrom", "start", "end", "exon"])
    except:
        raise Exception("Error! Could not read probes file. "
                        "Check file exists, is tab-delimited and has appropriate header?"
                        "Required fields: <chrom>, <start>, <end>, <exon_name>")

    probes["chrom"] = map(convert_chrom, probes.chrom)
    if probes["chrom"].dtype != np.int64:
        raise Exception("Error! Could not convert all probes to standard chromosome notation!")

    probes["probe_size"] = probes["end"] - probes["start"]

    if minsize:
        mask = probes.probe_size < minsize
        d = minsize - probes[mask]["probe_size"]
        left_pad = np.floor(d/2)
        right_pad = np.ceil(d/2)
        probes["start"][mask] = probes["start"][mask] - left_pad
        probes["end"][mask] = probes["end"][mask] + right_pad

    return probes

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("probes", help="probe file, bed format")
    parser.add_argument("input", help="Input file from frFAST mapper (*.hdf5)")
    parser.add_argument("output", help="Output hdf5 file containing RPKM values")
    parser.add_argument("--min-probe-size", default=10, type=int)
    parser.add_argument("--sampleID", help="Optional sample ID for file; default is to use input filename.", default=None)
    args = parser.parse_args()
    if not args.sampleID:
        args.sampleID = os.path.basename(args.input)

    t1 = time.time()
    try:
        infile = openFile(args.input, mode='r')
    except IOError:
        raise Exception("Error! Could not open %s for reading " % args.input)

    try:
        outfile = openFile(args.output, mode = 'w', title=str(args.sampleID))
    except IOError:
        raise Exception("Error! Could not obtain write handle for output file!")

    probes = load_probes(args.probes, minsize=args.min_probe_size)
    
    group = outfile.createGroup("/",'rpkm','rpkm')
    out_tbl = outfile.createTable(group,"rpkm",exon,"rpkm")
    
    # get total reads
    unique_reads_by_contig = infile.root.info_group.info_table.read(field='unique_reads')
    contig_names = infile.root.info_group.info_table.read(field='contig')

    total_reads = 0
    for contig,cnt in zip(contig_names,unique_reads_by_contig):
        if contig not in ["chr23","chr24"]:
            total_reads += int(cnt)
    
    read_t, calc_t, write_t = 0, 0, 0

    for chrom in xrange(1,25):

        chr_mask = probes.chrom == chrom
        exon_starts = probes[chr_mask]["start"].values
        exon_stops = probes[chr_mask]["end"].values
        exon_id = probes[chr_mask].index.values
        t2 = time.time()

        in_tbl = infile.root.read_starts._f_getChild("chr" + str(chrom))
        positions = in_tbl.read(field='pos')
        counts = in_tbl.read(field='read_starts')
        t3 = time.time()

        # RPKM calculation, very speed, such vectorized!
        a = np.searchsorted(positions, exon_starts + 1)
        b = np.searchsorted(positions, exon_stops + 1)
        ix = itertools.izip(a,b)
        read_cnt = np.array([np.sum(counts[j:k]) for j,k in ix])
        exon_bp = exon_stops-exon_starts
        rpkm = (10**9*(read_cnt)/(exon_bp))/(total_reads)
        t4 = time.time()

        out_data = np.empty([len(exon_id)],dtype='u4,u4,f8')
        out_data['f0'] = exon_id
        out_data['f1'] = read_cnt
        out_data['f2'] = rpkm
        out_tbl.append(out_data)
        
        t5 = time.time()
        read_t += t3 - t2
        calc_t += t4 - t3
        write_t += t5 - t4
    
    infile.close()
    outfile.close()
    print args.sampleID, "read: ", read_t, " calc: ", calc_t, " write: ", write_t, " total: ", time.time()-t1
