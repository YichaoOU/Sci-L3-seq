#!/usr/bin/env python

import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import os
import argparse
import pandas as pd
import subprocess
import glob
'''

In the working dir, run python step3.py {{jid}}



'''

def get_samples(jobID):

	out = subprocess.Popen('grep "Sample: " %s/*/*.out'%(jobID),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True).communicate()[0]
	out = out.split("\n")[:-1]

	df = pd.DataFrame([x.split() for x in out])
	df['sample'] = df[1]
	df['BC1'] = df[3].astype(int)
	df['BC2'] = df[5].astype(int)
	df['BC3'] = df[7].astype(int)
	for i in range(9):
		try:
			df = df.drop([i],axis=1)
		except:
			pass
	# print (df)
	df.index = df['sample']
	return df
	
def get_collision_rate(jobID):

	out = subprocess.Popen('grep "collision rate for sample" %s/*/*.out'%(jobID),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True).communicate()[0]
	out = out.split("\n")[:-1]

	df = pd.DataFrame([x.split() for x in out])
	df['sample'] = df[4]
	df['collision_rate'] = df[6].astype(float)
	for i in range(9):
		try:
			df = df.drop([i],axis=1)
		except:
			pass
	# print (df)
	df.index = df['sample']
	return df
	

def wccount(filename):
	out = subprocess.Popen("zcat %s | wc -l"%(filename),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True).communicate()[0]
	return int(out.partition(b' ')[0])/4
		

def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	mainParser.add_argument("-j","--jobID",  help="result_folder", required=True)
	mainParser.add_argument("-f","--input",  help="fastq.tsv", required=True)

	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def get_file(x):
	return glob.glob(x)[0]

def main():

	args = my_args()
	df = get_samples(args.jobID)
	df2 = get_collision_rate(args.jobID)
	df['collision_rate'] = df2['collision_rate']
	samples = pd.read_csv(args.input,sep="\t",header=None,index_col=2)
	# print samples
	samples = samples[0].to_dict()
	total_reads = []
	reads_left_step1 = []
	noRT =[]
	noBC3 = []
	noBC1_BC2 = []
	valid_reads = []
	for s in df['sample']:
		#total reads
		total_reads.append(wccount(samples[s]))
		# noRT
		noRT.append(wccount("%s/%s/%s.noRT.R1.fq.gz"%(args.jobID,s,s)))
		# noBC3
		noBC3.append(wccount("%s/%s/%s.noBC3.R1.fq.gz"%(args.jobID,s,s)))
		# reads_left_step1
		reads_left_step1.append(wccount("%s/%s/%s.R1.ordered.fastq.gz"%(args.jobID,s,s)))
		# noBC1_BC2
		noBC1_BC2.append(wccount("%s/%s/%s_barcode_demultiplexing/%s.junk.R1.fastq.gz"%(args.jobID,s,s,s)))
		# valid_reads
		valid_reads.append(wccount("%s/%s/%s_barcode_demultiplexing/%s.matched.R1.fastq.gz"%(args.jobID,s,s,s)))
	
	df['total_reads'] = total_reads
	df['reads_left_step1'] = reads_left_step1
	df['noRT'] = noRT
	df['noBC3'] = noBC3
	df['valid_reads'] = valid_reads
	df['noBC1_BC2'] = noBC1_BC2
	df['junk_reads'] = df['noRT'] +df['noBC3'] +df['noBC1_BC2']
	df['valid_reads_fraction'] = df['valid_reads']/df['total_reads']
	df['no_RT_fraction'] = df['noRT']/df['total_reads']
	df['no_BC3_fraction'] = df['noBC3']/df['total_reads']
	df['no_BC1BC2_fraction'] = df['noBC1_BC2']/df['total_reads']
	
	column_order = ['sample','collision_rate','valid_reads_fraction','total_reads','valid_reads','junk_reads','no_RT_fraction','noRT','no_BC3_fraction','noBC3','no_BC1BC2_fraction','noBC1_BC2','reads_left_step1','BC1','BC2','BC3']
	df[column_order].to_csv("%s/sample_QC.tsv"%(args.jobID),sep="\t",index=False)
	
if __name__ == "__main__":
	main()
	
		