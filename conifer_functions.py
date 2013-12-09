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
# CoNIFER is free software: you can redistribute it and/or modify
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

import csv
from tables import *
import numpy as np
import operator

class rpkm_value(IsDescription):
	probeID = UInt32Col(pos=0)
	rpkm = FloatCol(pos=1)

class probe(IsDescription):
	probeID = 	UInt32Col(pos=0)
	start = 	UInt32Col(pos=1) # start of probe
	stop = 		UInt32Col(pos=2) # stop of probe
	name = 		StringCol(20,pos=3)   # 20-character String

def chrInt2Str(chromosome_int):
	if int(chromosome_int) == 23:
		return 'chrX'
	elif int(chromosome_int) == 24:
		return 'chrY' 
	else:
		return 'chr' + str(chromosome_int)


def chrStr2Int(chromosome_str):
	chr = chromosome_str.replace('chr','')
	if chr == 'X':
		return 23
	elif chr == 'Y':
		return 24 
	else:
		return int(chr)

def parseLocString(locstr):
	try:
		chr,locstr = locstr.split(":")
		start, stop = locstr.split("-")
	except:
		chr, start, stop = locstr.split("\t")
	
	chr = chrStr2Int(chr)	
	start = int(start)
	stop = int(stop)
	return (chr,start,stop)	

def zrpkm(rpkm,median,sd):
	return (rpkm - median) / sd



class sample(IsDescription):
	sampleID = 	StringCol(100,pos=0) # 20-char string (sampleID)

def loadProbeList(CF_probe_filename):
	# Load data files
	probefile = open(CF_probe_filename, 'rb')
	s = csv.Sniffer()
	header = s.has_header(probefile.read(1024))
	probefile.seek(0)
	dialect = s.sniff(probefile.read(1024))
	probefile.seek(0)
	if header:
		r = csv.DictReader(probefile, dialect=dialect)
	else:
		r = csv.DictReader(probefile, dialect=dialect, fieldnames=['chr','start','stop','name'])
	
	probes = []
	probeID = 1
	for row in r:
		probes.append({'probeID': probeID, 'chr':chrStr2Int(row['chr']),'start':int(row['start']),'stop':int(row['stop']), 'name':row['name']})
		probeID +=1
	
	if len(probes) == 0:
		raise Exception("No probes in probe file")
		
	return probes


def loadRPKM(rpkm_filename):
	# test if rpkm_filename points to a HDF5 file or a txt file:
	if isHDF5File(rpkm_filename):
		rpkm_h5 = openFile(rpkm_filename, mode='r')
		return rpkm_h5.root.rpkm.rpkm.read(field="rpkm")
		rpkm_h5.close()
	else:
		return np.loadtxt(samples[s], dtype=np.float, delimiter="\t", skiprows=0, usecols=[2])
		
def export_sample(h5file_in,sample,probes,outfile_f):
	dt = np.dtype([('chr','|S10'),('start', '<u4'), ('stop', '<u4'), ('name', '|S20'),('SVDZRPKM',np.float)])
	for chr in h5file_in.root:
		if chr._v_title in ('probes','samples'):
			continue
		
		out_data = np.empty(len(probes[chr._v_title]),dtype=dt)
		out_data["SVDZRPKM"] = chr._f_getChild("sample_" + sample).read(field='rpkm')
		out_data["chr"] = np.repeat(chr._v_title,len(out_data))
		out_data["start"] = probes[chr._v_title]["start"]
		out_data["stop"] = probes[chr._v_title]["stop"]
		out_data["name"] = probes[chr._v_title]["name"]
		np.savetxt(outfile_f, out_data,fmt=["%s","%d","%d","%s","%f"], delimiter="\t")



def plotGenes(axis, rpkm_data, levels=5,x_pos=-2,text_pos='right',line_spacing=0.1,text_offset=0.25,data_range=None):
	from matplotlib.lines import Line2D
	counter = 0
	prev_gene = ""
	if data_range is not None:
		exon_set = rpkm_data.exons[data_range]
	else:
		exon_set = rpkm_data.exons
	for gene in exon_set["name"]:
		if gene == prev_gene:
			continue
		elif gene == 'None':
			continue
		start = np.min(np.where(exon_set["name"] == gene)) 
		stop = np.max(np.where(exon_set["name"] == gene)) + 1
		_ = axis.add_line(Line2D([start-0.5,stop-0.5],[x_pos - (counter * line_spacing),x_pos - (counter * line_spacing)],color=(102/255.,33/255.,168/255.,0.6),linewidth=5,linestyle='-',alpha=0.5,solid_capstyle='butt'))
		_ = axis.text(stop+text_offset, x_pos - (counter * line_spacing), gene, ha='left',va='center',fontsize=6)
		counter +=1
		prev_gene = gene
		if counter > 5:
			counter = 0


def plotGenomicCoords(plt, rpkm_data,fontsize=10,rotation=0):
	import operator
	import locale
	exon_set = rpkm_data.exons
	genomic_coords = np.array(map(operator.itemgetter("start"),exon_set))
	
	ticks = range(0,len(exon_set),len(exon_set)/5)
	
	ticks[-1] -= 1 # the last tick is going to be off the chart, so we estimate it as the second to last genomic coord.
	labels = [locale.format("%d", genomic_coords[i], grouping=True) for i in ticks if i < len(genomic_coords)]
	if rotation != 0:
		ha = "right"
	else:
		ha = "center"
	_ = plt.xticks(ticks,labels,fontsize=fontsize,rotation=rotation,ha=ha)


def plotRawData(axis, rpkm_data, color='r',linewidth=0.7):
	zero_stack = np.zeros(len(rpkm_data))
	positions = np.repeat(np.arange(0,len(rpkm_data)),3)
	logr = np.vstack([zero_stack,rpkm_data.flatten(),zero_stack]).transpose().ravel()
	axis.plot(positions,logr,color=color,marker=None,linewidth=1)

def getbkpoints(mask):
	bkpoints = np.nonzero(np.logical_xor(mask[0:-1],mask[1:]))[0]+1
	if mask[0] == 1:
		bkpoints = np.hstack([0,bkpoints])
	if mask[-1] == 1:
		bkpoints = np.hstack([bkpoints,len(mask)])
	return bkpoints.reshape(len(bkpoints)/2,2)

def mergeCalls(calls):
	if len(calls) == 0:
		return []
	
	out_calls = []
	calls=np.array(calls)[np.argsort(np.array(map(operator.itemgetter("start"),calls),dtype=np.int))]
	pstart = calls[0]["start"]
	pstop = calls[0]["stop"]
	for d in calls:
		if d["start"] <= pstop:
			pstop = max(d["stop"],pstop)
		else:
			out_calls.append({"sampleID": d["sampleID"], "chromosome": d["chromosome"], "start":pstart, "stop":pstop, "state": d["state"]})
			pstart = d["start"]
			pstop = d["stop"]
	
	out_calls.append({"sampleID": d["sampleID"], "chromosome": d["chromosome"], "start":pstart, "stop":pstop, "state": d["state"]})
	return out_calls

class rpkm_data:
	def __init__(self):
		self.rpkm = None
		self.samples = None
		self.exons = None
		self.isGenotype = False
		self.calls = []
		self.refined_calls = []
		
	def smooth(self, window = 15, padded = False): #todo, fix the padding here
		if self.isGenotype:
			print "Warning: the data in this rpkm_data container are single genotype values. Smoothing will have no effect!"
			return self.rpkm
		
		if window > 0:
			weightings=np.blackman(window)
			weightings = weightings/weightings.sum()
			smoothed_data = np.array([])
			for row in self.rpkm.transpose():
				smoothed = np.convolve(row, weightings)[(window-1)/2:-((window-1)/2)]
				if len(smoothed_data) == 0:
					smoothed_data  = smoothed
				else:
					smoothed_data  = np.vstack([smoothed_data,smoothed])
			
			self.rpkm = smoothed_data.transpose()
			return self.rpkm
		else:
			return self.rpkm
	
	def getSample(self, sampleIDs):
		sample_array = np.array(self.samples)
		if isinstance(sampleIDs,list):
			mask = np.zeros(len(sample_array),dtype=np.bool)
			for sampleID in sampleIDs:
				mask = np.logical_or(mask, sample_array == str(sampleID))
			
			return self.rpkm[:,mask]
		else:		
			mask = sample_array == str(sampleID)
			return self.rpkm[:,mask]
	
	def getSamples(self, sampleIDs):
		return self.getSample(sampleIDs)
	
	@property
	def shape(self):
		if self.isGenotype:
			return [len(self.samples), 1]
		else:
			return [len(self.samples), len(self.exons)]


class rpkm_reader:
	def __init__(self, rpkm_fn=None):
		"""Initialize an rpkm_reader instance. Specify the location of the data file"""
		
		if rpkm_fn == None:
			print "Must specify RPKM HDF5 file!"
			return 0
		# set up file access
		self.h5file = openFile(rpkm_fn, mode='r')
		self.sample_table = self.h5file.root.samples.samples
		
	def __del__(self):
		self.h5file.close()
	
	def getExonValuesByExons(self, chromosome, start_exon, stop_exon, sampleList=None,genotype=False):
		
		probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
		#table_rows = probe_tbl.getWhereList('(start >= %d) & (stop <= %d)' % (start,stop))
		start_exon = max(start_exon,0)
		stop_exon = min(stop_exon, probe_tbl.nrows)
		#print start_exon, stop_exon
		table_rows = np.arange(start_exon,stop_exon,1)
		data_tbl  = self.h5file.root._f_getChild("chr" + str(chromosome))
		
		if sampleList == None:
			num_samples = data_tbl._v_nchildren
			samples = data_tbl	
		else:
			num_samples = len(sampleList)
			samples = [data_tbl._f_getChild("sample_" + s) for s in sampleList]
		
		data = np.empty([num_samples,len(table_rows)],dtype=np.float)
		
		out_sample_list = []
		cnt = 0
		for sample_tbl in samples:
			d = sample_tbl.readCoordinates(table_rows,field="rpkm")
			data[cnt,:] = d
			cnt +=1
			out_sample_list.append(sample_tbl.title)
		
		d = rpkm_data()
		if genotype: # return average #todo-- implement median and SD?
			d.rpkm = data.transpose().mean(axis=0)
			d.isGenotype = True
		else: #return all data points
			d.rpkm = data.transpose()
		d.samples = out_sample_list
		d.exons = probe_tbl.readCoordinates(table_rows)
		
		return d
	
	def getExonValuesByRegion(self, chromosome, start=None, stop=None, sampleList=None,genotype=False):
		probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
		if (start is not None) and (stop is not None):
			table_rows = probe_tbl.getWhereList('(start >= %d) & (stop <= %d)' % (start,stop))
		else:
			table_rows = probe_tbl.getWhereList('(start >= 0) & (stop <= 1000000000)')
		
		data_tbl  = self.h5file.root._f_getChild("chr" + str(chromosome))
		
		if sampleList == None:
			num_samples = data_tbl._v_nchildren
			samples = data_tbl	
		else:
			num_samples = len(sampleList)
			samples = [data_tbl._f_getChild("sample_" + s) for s in sampleList]
		
		data = np.empty([num_samples,len(table_rows)],dtype=np.float)
		
		out_sample_list = []
		cnt = 0
		for sample_tbl in samples:
			d = sample_tbl.readCoordinates(table_rows,field="rpkm")
			data[cnt,:] = d
			cnt +=1
			out_sample_list.append(sample_tbl.title)
		
		d = rpkm_data()
		if genotype: # return average #todo-- implement median and SD?
			d.rpkm = data.transpose().mean(axis=0)
			d.isGenotype = True
		else: #return all data points
			d.rpkm = data.transpose()
		d.samples = out_sample_list
		d.exons = probe_tbl.readCoordinates(table_rows)
		
		return d
	
	def getSampleList(self,cohort=None,sex=None,ethnicity=None,custom=None):
		"""Return a list of available samples in the current data file. Specifying no arguments will return all available samples"""
		
		readWhereStr = ""
		if custom != None:
			readWhereStr = custom
		else:
			if cohort != None:
				if isinstance(cohort,list):
					for c in cohort:
						readWhereStr += "(cohort=='%s') | " % c
					readWhereStr = readWhereStr.strip(" |")
					readWhereStr += " & "
				else:
					readWhereStr += "(cohort=='%s') " % cohort
			if sex != None:
				if sex not in ['M','F']:	
					sex = sex.upper()[0]
				readWhereStr += " (sex=='%s') &" % sex
			if ethnicity != None:
				readWhereStr += " (ethnicity=='%s') &" % ethnicity
			
			readWhereStr = readWhereStr.strip(" &") # remove leading or trailing characters
		if readWhereStr != "":
			#print readWhereStr
			sampleIDs = self.sample_table.readWhere(readWhereStr,field='sampleID')
		else:
			sampleIDs = self.sample_table.read(field='sampleID')
		
		return sampleIDs
	
	def getExonIDs(self, chromosome, start, stop):
		probe_tbl = self.h5file.root.probes._f_getChild("probes_chr" + str(chromosome))
		exons = probe_tbl.getWhereList('(start >= %d) & (stop <= %d)' % (start,stop))
		return exons