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

python ../../../../src/step4_calculate_collision_rate.py --table KOK676_S3.total_number_reads.tsv --human human/KOK676_S3.R1.bed --mouse mouse/KOK676_S3.R1.bed

'''

def get_samples(jobID):

	out = subprocess.Popen('grep "Sample: " run8/*/*.out',stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True).communicate()[0]
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
	return df
	

def wccount(filename):
	out = subprocess.Popen("zcat %s | wc -l"%(filename),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True).communicate()[0]
	return int(out.partition(b' ')[0])/4
		

def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	mainParser.add_argument("--table",  help="barcode count table", required=True)
	mainParser.add_argument("--human",  help="human bamtobed r1", required=True)
	mainParser.add_argument("--mouse",  help="mouse bamtobed r1", required=True)
	mainParser.add_argument("--threshold",  help="threshold for collision", default=0.9,type=float)

	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def get_file(x):
	return glob.glob(x)[0]
	
def parse_file(x,mydict):
	with open(x) as f:
		for line in f:
			temp = line.split()[3].split(",")
			label = temp[1]+","+temp[2]+","+temp[3]
			mydict[label]+=1
	return mydict
			

def main():

	args = my_args()
	
	barcode = pd.read_csv(args.table,sep="\t")
	barcode = barcode[barcode['total_reads']>0]
	sample = args.table.split(".")[0]
	barcode['name'] = barcode['barcode_3']+","+barcode['barcode_2']+","+barcode['barcode_1']
	barcode.index = barcode['name']
	barcode['init'] = 0
	
	total_human = barcode['init'].to_dict()
	total_mouse = barcode['init'].to_dict()
	total_human = parse_file(args.human,total_human)
	total_mouse = parse_file(args.mouse,total_mouse)
	df1 = pd.DataFrame.from_dict(total_human,orient='index')
	df2 = pd.DataFrame.from_dict(total_mouse,orient='index')
	barcode['human_count']  = df1[0]
	barcode['mouse_count']  = df2[0]
	
	barcode['human_percent'] =  barcode['human_count'] / barcode['total_reads']
	barcode['mouse_percent'] =  barcode['mouse_count'] / barcode['total_reads']
	barcode['max_mapping_rate'] =  barcode[['human_percent','mouse_percent']].max(axis=1)
	barcode['collision'] = barcode['max_mapping_rate']<args.threshold
	# print (barcode.head())
	print ("collision rate for sample %s is %s "%(sample,float(barcode['collision'].sum())/barcode.shape[0]))
	barcode.to_csv("mapping_rate.tsv",sep="\t")
	barcode[['human_count','mouse_count','total_reads','human_percent','mouse_percent']].to_csv("for_collision_plot.tsv",sep=" ",index=True,header=False)
	

	
if __name__ == "__main__":
	main()
	
		