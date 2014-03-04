#######################################################################
#######################################################################
# CoNIFER: Copy Number Inference From Exome Reads
# Developed by Niklas Krumm (C) 2012
# nkrumm@gmail.com
# 
# homepage: http://conifer.sf.net
# This program is described in:
# Krumm et al. 2012. Genome Research. doi:10.1101/gr.138115.112
#
# This file is part of CoNIFER.
# CoNIFER is free software:  you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################
#######################################################################

import argparse
import os, sys, copy
import glob
import csv
import conifer_functions as cf
import operator
from tables import *
import numpy  as np

def CF_analyze(args):
	# do path/file checks:
	try: 
		# read probes table
		probe_fn = str(args.probes[0])
		probes = cf.loadProbeList(probe_fn)
		num_probes = len(probes)
		print '[INIT] Successfully read in %d probes from %s' % (num_probes, probe_fn)
	except IOError as e: 
		print '[ERROR] Cannot read probes file: ', probe_fn
		sys.exit(0)
	
	try: 
		svd_outfile_fn = str(args.output)
		h5file_out = openFile(svd_outfile_fn, mode='w')
		probe_group = h5file_out.createGroup("/","probes","probes")
	except IOError as e: 
		print '[ERROR] Cannot open SVD output file for writing: ', svd_outfile_fn
		sys.exit(0)
	
	if args.write_svals != "":
		sval_f = open(args.write_svals,'w')
	
	if args.plot_scree != "":
		try:
			import matplotlib
			matplotlib.use('Agg')
			import matplotlib.pyplot as plt
			import pylab as P
			from matplotlib.lines import Line2D
			from matplotlib.patches import Rectangle
		except:
			print "[ERROR] One or more of the required modules for plotting cannot be loaded! Are matplotlib and pylab installed?"
			sys.exit(0)
		
		plt.gcf().clear()
		fig = plt.figure(figsize=(10,5))
		ax = fig.add_subplot(111)
	
	rpkm_dir = str(args.rpkm_dir[0])
	rpkm_files = glob.glob(rpkm_dir + "/*")
	if len(rpkm_files) == 0:
		print '[ERROR] Cannot find any files in RPKM directory (or directory path is incorrect): ', rpkm_dir
		sys.exit(0)
	elif len(rpkm_files) == 1:
		print '[ERROR] Found only 1 RPKM file (sample). CoNIFER requires multiple samples (8 or more) to run. Exiting.'
		sys.exit(0)
	elif len(rpkm_files) < 8:
		print '[WARNING] Only found %d samples... this is less than the recommended minimum, and CoNIFER may not analyze this dataset correctly!' % len(rpkm_files)
	elif len(rpkm_files) <= int(args.svd):
		print '[ERROR] The number of SVD values specified (%d) must be less than the number of samples (%d). Either add more samples to the analysis or reduce the --svd parameter! Exiting.' % (len(rpkm_files), int(args.svd))
		sys.exit(0)
	else:
		print '[INIT] Found %d RPKM files in %s' % (len(rpkm_files), rpkm_dir)
	
	# read in samples names and generate file list
	samples = {}
	for f in rpkm_files:
		s = '.'.join(f.split('/')[-1].split('.')[0:-1])
		print "[INIT] Mapping file to sampleID: %s --> %s" % (f, s)
		samples[s] = f
	
	#check uniqueness and total # of samples
	if len(set(samples)) != len(set(rpkm_files)):
		print '[ERROR] Could not successfully derive sample names from RPKM filenames. There are probably non-unique sample names! Please rename files using <sampleID>.txt format!'
		sys.exit(0)
	
	# LOAD RPKM DATA
	RPKM_data = np.zeros([num_probes,len(samples)], dtype=np.float)
	failed_samples = 0
	
	for i,s in enumerate(samples.keys()):
		t = cf.loadRPKM(samples[s])
		if len(t) != num_probes:
			print "[WARNING] Number of RPKM values for %s in file %s does not match number of defined probes in %s. **This sample will be dropped from analysis**!" % (s, samples[s], probe_fn)
			_ = samples.pop(s)
			failed_samples += 1
		else:
			RPKM_data[:,i] = t
			print "[INIT] Successfully read RPKM data for sampleID: %s" % s	
	
	RPKM_data = RPKM_data[:,0:len(samples)]
	print "[INIT] Finished reading RPKM files. Total number of samples in experiment: %d (%d failed to read properly)" % (len(samples), failed_samples)
	
	if len(samples) <= int(args.svd):
		print '[ERROR] The number of SVD values specified (%d) must be less than the number of samples (%d). Either add more samples to the analysis or reduce the --svd parameter! Exiting.' % (int(args.svd), len(samples))
		sys.exit(0)
	
	# BEGIN 
	chrs_to_process = set(map(operator.itemgetter("chr"),probes))
	chrs_to_process_str = ', '.join([cf.chrInt2Str(c) for c in chrs_to_process])
	print '[INIT] Attempting to process chromosomes: ', chrs_to_process_str
	
	if args.all_chromosomes:
		print "[RUNNING: all chromosomes] Found %d probes" % len(probes)
		print "[RUNNING: all chromosomes] Calculating median RPKM"
		median = np.median(RPKM_data,1)
		sd = np.std(RPKM_data,1)
		probe_mask = median >= float(args.min_rpkm)
		print "[RUNNING: all chromosomes] Masking %d probes with median RPKM < %f" % (np.sum(probe_mask==False), float(args.min_rpkm))
		RPKM_data = RPKM_data[probe_mask, :]
		num_chr_probes = np.sum(probe_mask)
		print "[RUNNING: all chromosomes] Calculating ZRPKM scores..."
		RPKM_data = np.apply_along_axis(cf.zrpkm, 0, RPKM_data[probe_mask], median[probe_mask], sd[probe_mask])
		print "[RUNNING: all chromosomes] SVD decomposition..."
		components_removed = int(args.svd)
		
		U, S, Vt = np.linalg.svd(RPKM_data,full_matrices=False)
		new_S = np.diag(np.hstack([np.zeros([components_removed]),S[components_removed:]]))
		
		if args.write_svals != "":
			sval_f.write('chr' + str(chr) + '\t' + '\t'.join([str(_i) for _i in S]) + "\n")
		
		if args.plot_scree != "":
			ax.plot(S, label='chr' + str(chr),lw=0.5)
		
		# reconstruct data matrix
		RPKM_data = np.dot(U, np.dot(new_S, Vt))
	
	for chr in chrs_to_process:
		print "[RUNNING: chr%d] Now on: %s" %(chr, cf.chrInt2Str(chr))
		chr_group_name = "chr%d" % chr
		chr_group = h5file_out.createGroup("/",chr_group_name,chr_group_name)
		
		chr_mask = np.array(map(operator.itemgetter("chr"),probes)) == chr
		
		if not args.all_chromosomes:
			chr_probes = filter(lambda i: i["chr"] == chr, probes)
			num_chr_probes = len(chr_probes)
			start_probeID = chr_probes[0]['probeID']
			stop_probeID = chr_probes[-1]['probeID']
			print "[RUNNING: chr%d] Found %d probes; probeID range is [%d-%d]" % (chr, len(chr_probes), start_probeID-1, stop_probeID) # probeID is 1-based and slicing is 0-based, hence the start_probeID-1 term
			
			rpkm = RPKM_data[start_probeID:stop_probeID,:]
			
			print "[RUNNING: chr%d] Calculating median RPKM" % chr
			median = np.median(rpkm,1)
			sd = np.std(rpkm,1)
			probe_mask = median >= float(args.min_rpkm)
			print "[RUNNING: chr%d] Masking %d probes with median RPKM < %f" % (chr, np.sum(probe_mask==False), float(args.min_rpkm))
			
			rpkm = rpkm[probe_mask, :]
			num_chr_probes = np.sum(probe_mask)
			
			if num_chr_probes <= len(samples):
				print "[ERROR] This chromosome has fewer informative probes than there are samples in the analysis! There are probably no mappings on this chromosome. Please remove these probes from the probes.txt file"
				sys.exit(0)

			print "[RUNNING: chr%d] Calculating ZRPKM scores..." % chr
			rpkm = np.apply_along_axis(cf.zrpkm, 0, rpkm, median[probe_mask], sd[probe_mask])
			
			# svd transform
			print "[RUNNING: chr%d] SVD decomposition..." % chr
			components_removed = int(args.svd)
			
			U, S, Vt = np.linalg.svd(rpkm,full_matrices=False)
			new_S = np.diag(np.hstack([np.zeros([components_removed]),S[components_removed:]]))
			
			if args.write_svals != "":
				sval_f.write('chr' + str(chr) + '\t' + '\t'.join([str(_i) for _i in S]) + "\n")
			
			if args.plot_scree != "":
				ax.plot(S, label='chr' + str(chr),lw=0.5)
			
			# reconstruct data matrix
			rpkm = np.dot(U, np.dot(new_S, Vt))

			probeIDs = np.array(map(operator.itemgetter("probeID"),chr_probes))[probe_mask]
			probe_starts = np.array(map(operator.itemgetter("start"),chr_probes))[probe_mask]
			probe_stops = np.array(map(operator.itemgetter("stop"),chr_probes))[probe_mask]	
			gene_names =  np.array(map(operator.itemgetter("name"),chr_probes))[probe_mask]	
		else:
			num_chr_probes = np.sum(chr_mask & probe_mask)
			probeIDs = np.array(map(operator.itemgetter("probeID"),probes))[chr_mask & probe_mask]
			probe_starts = np.array(map(operator.itemgetter("start"),probes))[chr_mask & probe_mask]
			probe_stops = np.array(map(operator.itemgetter("stop"),probes))[chr_mask & probe_mask]	
			gene_names =  np.array(map(operator.itemgetter("name"),probes))[chr_mask & probe_mask]	

			rpkm = RPKM_data[chr_mask & probe_mask, :]
		
		dt = np.dtype([('probeID',np.uint32),('start',np.uint32),('stop',np.uint32), ('name', np.str_, 20)])
		
		out_probes = np.empty(num_chr_probes,dtype=dt)
		out_probes['probeID'] = probeIDs
		out_probes['start'] = probe_starts
		out_probes['stop'] = probe_stops
		out_probes['name'] = gene_names
		probe_table = h5file_out.createTable(probe_group,"probes_chr%d" % chr,cf.probe,"chr%d" % chr)
		probe_table.append(out_probes)
		
		
		# save to HDF5 file
		print "[RUNNING: chr%d] Saving SVD-ZRPKM values" % chr
		
		for i,s in enumerate(samples):
			out_data = np.empty(num_chr_probes,dtype='u4,f8')
			out_data['f0'] = probeIDs
			out_data['f1'] = rpkm[:,i]
			sample_tbl = h5file_out.createTable(chr_group,"sample_" + str(s),cf.rpkm_value,"%s" % str(s))
			sample_tbl.append(out_data)
	
	
	print "[RUNNING] Saving sampleIDs to file..."
	sample_group = h5file_out.createGroup("/","samples","samples")
	sample_table = h5file_out.createTable(sample_group,"samples",cf.sample,"samples")
	dt = np.dtype([('sampleID',np.str_,100)])
	out_samples = np.empty(len(samples.keys()),dtype=dt)
	out_samples['sampleID'] = np.array(samples.keys())
	sample_table.append(out_samples)
	
	
	if args.write_sd != "":
		print "[RUNNING] Calculating standard deviations for all samples (this can take a while)..."
		
		sd_file = open(args.write_sd,'w')
		
		for i,s in enumerate(samples):
			# collect all SVD-ZRPKM values
			count = 1
			for chr in chrs_to_process:
				if count == 1:
					sd_out = h5file_out.root._f_getChild("chr%d" % chr)._f_getChild("sample_%s" % s).read(field="rpkm").flatten()
				else:
					sd_out = np.hstack([sd_out,out.h5file_out.root._f_getChild("chr%d" % chr)._f_getChild("sample_%s" % s).read(field="rpkm").flatten()])
				
				sd = np.std(sd_out)
			sd_file.write("%s\t%f\n" % (s,sd))
		
		sd_file.close()
	
	if args.plot_scree != "":
		plt.title("Scree plot")
		if len(samples) < 50:
			plt.xlim([0,len(samples)])
			plt.xlabel("S values")
		else:
			plt.xlim([0,50])
			plt.xlabel("S values (only first 50 plotted)")
		plt.ylabel("Relative contributed variance")		
		plt.savefig(args.plot_scree)
	
	print "[FINISHED]"
	h5file_out.close()
	sys.exit(0)

def CF_export(args):
	try: 
		h5file_in_fn = str(args.input)
		h5file_in = openFile(h5file_in_fn, mode='r')
	except IOError as e: 
		print '[ERROR] Cannot open CoNIFER input file for reading: ', h5file_in_fn
		sys.exit(0)	
	
	# read probes
	probes = {}
	for probes_chr in h5file_in.root.probes:
		probes[probes_chr.title] = probes_chr.read()
	
	if args.sample =='all':
		all_samples = list(h5file_in.root.samples.samples.read(field="sampleID"))
		
		out_path = os.path.abspath(args.output)
		
		print "[INIT] Preparing to export all samples (%d samples) to %s" % (len(all_samples), out_path)
		for sample in all_samples:
			try:
				outfile_fn = out_path + "/" + sample + ".bed"
				outfile_f = open(outfile_fn,'w')
			except IOError as e:
				print '[ERROR] Cannot open output file for writing: ', outfile_fn
				sys.exit(0)
			print "[RUNNING] Exporting %s" % sample
			
			cf.export_sample(h5file_in,sample,probes,outfile_f)
			outfile_f.close()
	
	elif len(args.sample) == 1:
		out_path = os.path.abspath(args.output)
		sample = args.sample[0]
		print "[INIT] Preparing to export sampleID %s to %s" % (args.sample[0], out_path)
		try:
			if os.path.isdir(out_path):
				outfile_fn = out_path + "/" + sample + ".bed"
			else:
				outfile_fn = out_path
			outfile_f = open(outfile_fn,'w')
		except IOError as e:
			print '[ERROR] Cannot open output file for writing: ', outfile_fn
			sys.exit(0)
		print "[RUNNING] Exporting %s to %s" % (sample, outfile_fn)
		
		cf.export_sample(h5file_in,sample,probes,outfile_f)
		outfile_f.close()
	
	else:
		out_path = os.path.abspath(args.output)
		print "[INIT] Preparing to export %d samples to %s" % (len(args.sample), out_path)
		for sample in args.sample:
			try:
				if os.path.isdir(out_path):
					outfile_fn = out_path + "/" + sample + ".bed"
				else:
					outfile_fn = out_path
				outfile_f = open(outfile_fn,'w')
			except IOError as e:
				print '[ERROR] Cannot open output file for writing: ', outfile_fn
				sys.exit(0)
			print "[RUNNING] Exporting %s to %s" % (sample, outfile_fn)
			
			cf.export_sample(h5file_in,sample,probes,outfile_f)
			outfile_f.close()		
	sys.exit(0)

def CF_call(args):
	try: 
		h5file_in_fn = str(args.input)
		h5file_in = openFile(h5file_in_fn, mode='r')
	except IOError as e: 
		print '[ERROR] Cannot open CoNIFER input file for reading: ', h5file_in_fn
		sys.exit(0)		
	
	try: 
		callfile_fn = str(args.output)
		callfile_f = open(callfile_fn, mode='w')
	except IOError as e: 
		print '[ERROR] Cannot open output file for writing: ', callfile_fn
		sys.exit(0)
	
	chrs_to_process = []
	for chr in h5file_in.root:
		if chr._v_title not in ('probes','samples'):
			chrs_to_process.append(chr._v_title.replace("chr",""))
	
	h5file_in.close()
	
	print '[INIT] Initializing caller at threshold = %f' % (args.threshold)
	
	r = cf.rpkm_reader(h5file_in_fn)
	
	all_calls = []
	
	for chr in chrs_to_process:
		print '[RUNNING] Now processing chr%s' % chr
		data = r.getExonValuesByRegion(chr)
		
		#raw_data = copy.copy(data)
		_ = data.smooth()
		
		mean= np.mean(data.rpkm,axis=1)
		sd =  np.std(data.rpkm,axis=1)
		
		for sample in r.getSampleList():
			sample_data = data.getSample([sample]).flatten()
			#sample_raw_data = raw_data.getSample([sample]).flatten()
			
			dup_mask = sample_data >= args.threshold
			del_mask = sample_data <= -1*args.threshold
			
			dup_bkpoints = cf.getbkpoints(dup_mask) #returns exon coordinates for this chromosome (numpy array coords)
			del_bkpoints = cf.getbkpoints(del_mask)
			
			
			dups = []
			for start,stop in dup_bkpoints:
				try: new_start =  np.max(np.where(sample_data[:start] < (mean[:start] + 3*sd[:start])))
				except ValueError: new_start = 0
				try: new_stop = stop + np.min(np.where(sample_data[stop:] < (mean[stop:] + 3*sd[stop:])))
				except ValueError: new_stop = data.shape[1]-1
				dups.append({"sampleID":sample,"chromosome":  cf.chrInt2Str(chr), "start":data.exons[new_start]["start"], "stop": data.exons[new_stop]["stop"], "state": "dup"})
			
			dels = []
			for start,stop in del_bkpoints:	
				try: new_start =  np.max(np.where(sample_data[:start] > (-1*mean[:start] - 3*sd[:start])))
				except ValueError: new_start = 0
				try: new_stop = stop + np.min(np.where(sample_data[stop:] > (-1*mean[stop:] - 3*sd[stop:])))
				except ValueError: new_stop = data.shape[1]-1
				dels.append({"sampleID":sample,"chromosome": cf.chrInt2Str(chr), "start":data.exons[new_start]["start"], "stop": data.exons[new_stop]["stop"], "state": "del"})
			
			dels = cf.mergeCalls(dels) #merges overlapping calls
			dups = cf.mergeCalls(dups)
			
			#print sampleID, len(dels), len(dups)
			
			all_calls.extend(list(dels))
			all_calls.extend(list(dups))
	
	# print calls to file
	header = ['sampleID','chromosome','start','stop','state']
	
	callfile_f.write('\t'.join(header) + "\n")
	for call in all_calls:
		print "%s\t%s\t%d\t%d\t%s" % (call["sampleID"], call["chromosome"], call["start"], call["stop"], call["state"])
		callfile_f.write("%s\t%s\t%d\t%d\t%s\n" % (call["sampleID"], call["chromosome"], call["start"], call["stop"], call["state"]))
	
	sys.exit(0)

def CF_plot(args):
	try:
		import locale
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		import pylab as P
		from matplotlib.lines import Line2D
		from matplotlib.patches import Rectangle
		_ = locale.setlocale(locale.LC_ALL, '')
	except:
		print "[ERROR] One or more of the required modules for plotting cannot be loaded! Are matplotlib and pylab installed?"
		sys.exit(0)
		
		
	chr, start, stop = cf.parseLocString(args.region)
	
	r = cf.rpkm_reader(args.input)
	
	data = r.getExonValuesByRegion(chr,start,stop)
	_ = data.smooth()
	
	plt.gcf().clear()
	fig = plt.figure(figsize=(10,5))
	ax = fig.add_subplot(111)
	
	
	ax.plot(data.rpkm, linewidth = 0.3, c='k')
	
	
	if args.sample != 'none':
		cnt = 1
		coloriter = iter(['r','b','g','y'])
		for sample in args.sample:
			try:
				color, sampleID = sample.split(":")
			except:
				color =coloriter.next()
				sampleID = sample
			
			ax.plot(data.getSample([sampleID]), linewidth = 1, c=color, label = sampleID)
			
			if cnt == 1:
				cf.plotRawData(ax, r.getExonValuesByRegion(chr,start,stop,sampleList=[sampleID]).getSample([sampleID]),color=color)
			cnt +=1
		plt.legend(prop={'size':10},frameon=False)
		
	cf.plotGenes(ax, data)
	cf.plotGenomicCoords(plt,data)
	plt.xlim(0,data.shape[1])
	plt.ylim(-3,3)
	
	plt.title("%s: %s - %s" % (cf.chrInt2Str(chr),locale.format("%d",start, grouping=True),locale.format("%d",stop, grouping=True)))
	plt.xlabel("Position")
	plt.ylabel("SVD-ZRPKM Values")
	
	plt.savefig(args.output)
	
	sys.exit(0)

def CF_plotcalls(args):
	try:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		import pylab as P
		from matplotlib.lines import Line2D
		from matplotlib.patches import Rectangle
	except:
		print "[ERROR] One or more of the required modules for plotting cannot be loaded! Are matplotlib and pylab installed?"
		sys.exit(0)
	
	import locale	
	try:
		_ = locale.setlocale(locale.LC_ALL, 'en_US')
	except:
		_ = locale.setlocale(locale.LC_ALL, '')
	
	try: 
		callfile_fn = str(args.calls)
		callfile_f = open(callfile_fn, mode='r')
	except IOError as e: 
		print '[ERROR] Cannot open call file for reading: ', callfile_fn
		sys.exit(0)
	
	all_calls = []
	header = callfile_f.readline()
	
	for line in callfile_f:
		sampleID, chr, start, stop, state = line.strip().split()
		chr = cf.chrStr2Int(chr)
		all_calls.append({"chromosome":int(chr), "start":int(start), "stop":int(stop), "sampleID":sampleID})
	
	r = cf.rpkm_reader(args.input)
	
	for call in all_calls:
		chr = call["chromosome"]
		start = call["start"]
		stop = call["stop"]
		sampleID = call["sampleID"]
		
		exons = r.getExonIDs(chr,int(start),int(stop))
		
		
		window_start = max(exons[0]-args.window,0)
		window_stop = exons[-1]+args.window
		
		data = r.getExonValuesByExons(chr,window_start, window_stop)
		_ = data.smooth()
		
		plt.gcf().clear()
		fig = plt.figure(figsize=(10,5))
		ax = fig.add_subplot(111)
		
		
		ax.plot(data.rpkm, linewidth = 0.3, c='k')
		
		
		ax.plot(data.getSample([sampleID]), linewidth = 1, c='r', label = sampleID)
		cf.plotRawData(ax, r.getExonValuesByExons(chr,window_start, window_stop,sampleList=[sampleID]).getSample([sampleID]),color='r')
		
		plt.legend(prop={'size':10},frameon=False)
		
		cf.plotGenes(ax, data)
		cf.plotGenomicCoords(plt,data)
		
		exon_start = np.where(data.exons["start"] == start)[0][0]
		exon_stop = np.where(data.exons["stop"] == stop)[0][0]
		_ = ax.add_line(matplotlib.lines.Line2D([exon_start,exon_stop],[2,2],color='k',lw=6,linestyle='-',alpha=1,solid_capstyle='butt'))
		
		_ = plt.xlim(0,data.shape[1])
		_ = plt.ylim(-3,3)
		
		plt.title("%s: %s - %s" % (cf.chrInt2Str(chr),locale.format("%d",start, grouping=True),locale.format("%d",stop, grouping=True)))
		plt.xlabel("Position")
		plt.ylabel("SVD-ZRPKM Values")
		outfile = "%s_%d_%d_%s.png" % (cf.chrInt2Str(chr), start, stop, sampleID)
		plt.savefig(args.outputdir + "/" + outfile)

def CF_bam2RPKM(args):
	try:
		import pysam
	except:
		print '[ERROR] Cannot load pysam module! Make sure it is insalled'
		sys.exit(0)
	try: 
		# read probes table
		probe_fn = str(args.probes[0])
		probes = cf.loadProbeList(probe_fn)
		num_probes = len(probes)
		print '[INIT] Successfully read in %d probes from %s' % (num_probes, probe_fn)
	except IOError as e: 
		print '[ERROR] Cannot read probes file: ', probe_fn
		sys.exit(0)
	
	try:
		rpkm_f = open(args.output[0],'w')
	except IOError as e:
		print '[ERROR] Cannot open rpkm file for writing: ', args.output
		sys.exit(0)
	
	print "[RUNNING] Counting total number of reads in bam file..."
	total_reads = float(pysam.view("-c", args.input[0])[0].strip("\n"))
	print "[RUNNING] Found %d reads" % total_reads
	
	f = pysam.Samfile(args.input[0], "rb" )	
	
	if not f._hasIndex():
		print "[ERROR] No index found for bam file (%s)!\n[ERROR] You must first index the bam file and include the .bai file in the same directory as the bam file!" % args.input[0]
		sys.exit(0)
    
	# will be storing values in these arrays
	readcount = np.zeros(num_probes)
	exon_bp = np.zeros(num_probes)
	probeIDs = np.zeros(num_probes)
	counter = 0
	
	# detect contig naming scheme here # TODO, add an optional "contigs.txt" file or automatically handle contig naming
	bam_contigs = f.references
	probes_contigs = [str(p) for p in set(map(operator.itemgetter("chr"),probes))]
	
	probes2contigmap = {}
	
	for probes_contig in probes_contigs:
		if probes_contig in bam_contigs:
			probes2contigmap[probes_contig] = probes_contig
		elif cf.chrInt2Str(probes_contig) in bam_contigs:
			probes2contigmap[probes_contig] = cf.chrInt2Str(probes_contig)
		elif cf.chrInt2Str(probes_contig).replace("chr","") in bam_contigs:
			probes2contigmap[probes_contig] = cf.chrInt2Str(probes_contig).replace("chr","")
		else:
			print "[ERROR] Could not find contig '%s' from %s in bam file! \n[ERROR] Perhaps the contig names for the probes are incompatible with the bam file ('chr1' vs. '1'), or unsupported contig naming is used?" % (probes_contig, probe_fn)
			sys.exit(0)
	
	print "[RUNNING] Calculating RPKM values..."
	
	# loop through each probe	
	for p in probes:
		
		# f.fetch is a pysam method and returns an iterator for reads overlapping interval
		
		p_chr = probes2contigmap[str(p["chr"])]
		
		p_start = p["start"]
		p_stop = p["stop"]
		try:
			iter = f.fetch(p_chr,p_start,p_stop)
		except:
			print "[ERROR] Could not retrieve mappings for region %s:%d-%d. Check that contigs are named correctly and the bam file is properly indexed" % (p_chr,p_start,p_stop)
			sys.exit(0)
		
		for i in iter:
			if i.pos+1 >= p_start: #this checks to make sure a read actually starts in an interval
				readcount[counter] += 1
		
		exon_bp[counter] = p_stop-p_start
		probeIDs[counter] = counter +1 #probeIDs are 1-based
		counter +=1
	
	#calcualte RPKM values for all probes
	rpkm = (10**9*(readcount)/(exon_bp))/(total_reads)
	
	out = np.vstack([probeIDs,readcount,rpkm])
	
	np.savetxt(rpkm_f,out.transpose(),delimiter='\t',fmt=['%d','%d','%f'])
	
	rpkm_f.close()

	

VERSION = "targetedseq-dev-0.1"
parser = argparse.ArgumentParser(prog="CoNIFER", description="This is CoNIFER %s (Copy Number Inference From Exome Reads), designed to detect and genotype CNVs and CNPs from exome sequence read-depth data. See Krumm et al., Genome Research (2012) doi:10.1101/gr.138115.112 \nNiklas Krumm, 2012\n nkrumm@uw.edu" % VERSION)
parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
subparsers = parser.add_subparsers(help='Command to be run.')

# rpkm command
rpkm_parser= subparsers.add_parser('rpkm', help='Create an RPKM file from a BAM file')
rpkm_parser.add_argument('--probes',action='store', required=True, metavar='/path/to/probes_file.txt',  nargs=1,help="Probe definition file")
rpkm_parser.add_argument('--input',action='store', required=True, metavar='sample.bam',nargs=1,  help="Aligned BAM file")
rpkm_parser.add_argument('--output',action='store', required=True, metavar='sample.rpkm.txt',nargs=1,  help="RPKM file to write")
rpkm_parser.set_defaults(func=CF_bam2RPKM)

# analyze command
analyze_parser= subparsers.add_parser('analyze', help='Basic CoNIFER analysis. Reads a directory of RPKM files and a probe list and outputs a HDF5 file containing SVD-ZRPKM values.')
analyze_parser.add_argument('--probes',action='store', required=True, metavar='/path/to/probes_file.txt',  nargs=1,help="Probe definition file")
analyze_parser.add_argument('--rpkm_dir',action='store', required=True, metavar='/path/to/folder/containing/rpkm_files/',nargs=1,  help="Location of folder containing RPKM files. Folder should contain ONLY contain RPKM files, and all readable RPKM files will used in analysis.")
analyze_parser.add_argument('--output','-o', required=True, metavar='/path/to/output_file.hdf5', type=str, help="Output location of HDF5 file to contain SVD-ZRPKM values")
analyze_parser.add_argument('--svd', metavar='12', type=int, nargs='?', default = 12,help="Number of components to remove")
analyze_parser.add_argument('--min_rpkm', metavar='1.00', type=float, nargs="?", default = 1.00,help="Minimum population median RPKM per probe.")
analyze_parser.add_argument('--write_svals', metavar='SingularValues.txt', type=str, default= "", help="Optional file to write singular values (S-values). Used for Scree Plot.")
analyze_parser.add_argument('--plot_scree', metavar='ScreePlot.png', type=str, default= "", help="Optional graphical scree plot. Requires matplotlib.")
analyze_parser.add_argument('--write_sd', metavar='StandardDeviations.txt', type=str, default= "", help="Optional file with sample SVD-ZRPKM standard deviations. Used to filter noisy samples.")
analyze_parser.add_argument('--all_chromosomes', action="store_true", default=False, help="Run all chromosomes in one SVD transformation. Useful for targeted sequencing panels.")
analyze_parser.set_defaults(func=CF_analyze)

# export command
export_parser= subparsers.add_parser('export', help='Export SVD-ZRPKM values from a HDF5 file to bed or vcf format.')
export_parser.add_argument('--input','-i',action='store', required=True, metavar='CoNIFER_SVD.hdf5',help="HDF5 file from CoNIFER 'analyze' step")
export_parser.add_argument('--output','-o',action='store', required=False, default='.', metavar='output.bed',help="Location of output file[s]. When exporting multiple samples, top-level directory of this path will be used.")
export_parser.add_argument('--sample','-s',action='store', required=False, metavar='sampleID', default='all', nargs="+",help="SampleID or comma-separated list of sampleIDs to export. Default: export all samples")
export_parser.set_defaults(func=CF_export)

# plot command
plot_parser= subparsers.add_parser('plot', help='Plot SVD-ZRPKM values using matplotlib')
plot_parser.add_argument('--input','-i',action='store', required=True, metavar='CoNIFER_SVD.hdf5',help="HDF5 file from CoNIFER 'analyze' step")
plot_parser.add_argument('--region',action='store', required=True, metavar='chr#:start-stop',help="Region to plot")
plot_parser.add_argument('--output',action='store', required=True, metavar='image.png',help="Output path and filetype. PDF, PNG, PS, EPS, and SVG are supported.")
plot_parser.add_argument('--sample',action='store', required=False, metavar='sampleID',nargs="+",default='none',help="Sample[s] to highlight. The following optional color spec can be used: <color>:<sampleID>. Available colors are r,b,g,y,c,m,y,k. The unsmoothed SVD-ZRPKM values for the first sample in this list will be drawn. Default: No samples highlighted.")
plot_parser.set_defaults(func=CF_plot)

# make calls command
call_parser= subparsers.add_parser('call', help='Very rudimentary caller for CNVs using SVD-ZRPKM thresholding.')
call_parser.add_argument('--input','-i',action='store', required=True, metavar='CoNIFER_SVD.hdf5',help="HDF5 file from CoNIFER 'analyze' step")
call_parser.add_argument('--output',action='store', required=True, metavar='calls.txt',help="Output file for calls")
call_parser.add_argument('--threshold', metavar='1.50', type=float, nargs='?', required=False, default = 1.50,help="+/- Threshold for calling (minimum SVD-ZRPKM)")
call_parser.set_defaults(func=CF_call)

# plotcalls command
plotcalls_parser= subparsers.add_parser('plotcalls', help='Make basic plots from call file from "call" command.')
plotcalls_parser.add_argument('--input','-i',action='store', required=True, metavar='CoNIFER_SVD.hdf5',help="HDF5 file from CoNIFER 'analyze' step")
plotcalls_parser.add_argument('--calls',action='store', required=True, metavar='calls.txt',help="File with calls from 'call' command.")
plotcalls_parser.add_argument('--outputdir',action='store', required=True, metavar='/path/to/directory',help="Output directory for plots")
plotcalls_parser.add_argument('--window',action='store', required=False, metavar='50',default=50,help="In exons, the amount of padding to plot around each call")
plotcalls_parser.set_defaults(func=CF_plotcalls)

args = parser.parse_args()
args.func(args)
