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
import itertools
'''
	this script accepts and opens a R1 file, a R2 file, an RT primer sequence, and a SSS barcode list; then reorders the read pairs, and extracts and matches the nearest sss barcode and attach the sss barcode to R1 and R2 read names.
	the output files include ordered R1 and R2 files, noRTprimer.R1.fq.gz and noRTprimer.R2.fq.gz
	
python ../../../src/step2_demultiplexing.py -r1 KOK675_S2.noRTprimer.R1.fq.gz -r2 KOK675_S2.noRTprimer.R2.fq.gz --sample_ID test --barcode_1_list barcode_1.list --barcode_2_list barcode_2.list --barcode_3_list SSS_barcode.list --UMI_length 4 --bc3_length 6 --sp2_length 20 --bc2_length 7 --sp1_length 6 --bc1_length 8
	
'''

def dict3d_to_df(x):
	bc3 = []
	bc2 = []
	bc1 = []
	counts = []
	for i in x:
		for j in x[i]:
			for k in x[i][j]:
				bc3.append(i)
				bc2.append(j)
				bc1.append(k)
				counts.append(x[i][j][k])
	df = pd.DataFrame()
	df['barcode_3'] = bc3
	df['barcode_2'] = bc2
	df['barcode_1'] = bc1
	df['total_reads'] = counts
	return df
	

def k_mer_distance(k,error):
	out_dict={}
	bases=['A','T','G','C']
	k_mer = [''.join(p) for p in itertools.product(bases, repeat=k)]
	for ii in range(len(k_mer)):
		i = k_mer[ii]
		out_dict[i]={}
		for jj in range(ii,len(k_mer)):
			j = k_mer[jj]
			dist = distance(i, j)
			if dist <= error:
				out_dict[i][j] = distance
				out_dict[j][i] = distance
	return out_dict

def get_barcode_dict(bc1,bc2,bc3):
	out = {}
	for i in bc3:
		out[i]={}
		for j in bc2:
			out[i][j]={}
			for k in bc1:
				out[i][j][k]=0
	return out

def find_dist(kmer_dict,target_barcode,barcode_list):
	for x in barcode_list:
		try:
			kmer_dict[x][target_barcode] 
		except:
			return False
	return True


def output_to_fastq_gz(file_name,list_of_lines):
	f = gzip.open(file_name,"wb")
	[f.write(x) for x in list_of_lines]
	f.close()

def sci_l3_demultiplexing(Read1,Read2,label, barcode_1_list, barcode_2_list, barcode_3_list, BC1_error, BC2_error, BC3_error,UMI_length,bc3_length,sp2_length,bc2_length,sp1_length,bc1_length):
	"""general utils for demultiplexing umi-bc3-sp2-bc2-sp1-bc1-xxxx R1 read in PE mode
	
	
	"""
	
	output_folder = label+"_barcode_demultiplexing"
	os.system("mkdir -p %s"%(output_folder))
	barcode_dict = get_barcode_dict(barcode_1_list,barcode_2_list,barcode_3_list) # used to count reads
	# barcode_dict_R2 = get_barcode_dict(barcode_1_list,barcode_2_list,barcode_3_list)
	junk_list_R1 = []
	junk_list_R2 = []
	matched_list_R1 = []
	matched_list_R2 = []
	f1 = gzip.open(Read1)
	f2 = gzip.open(Read2)	
	
	bc3_kmer = k_mer_distance(bc3_length,BC3_error)
	bc2_kmer = k_mer_distance(bc2_length,BC2_error)
	bc1_kmer = k_mer_distance(bc1_length,BC1_error)

	bc3_start = UMI_length
	bc3_end = UMI_length+bc3_length
	bc2_start = UMI_length+bc3_length+sp2_length
	bc2_end = UMI_length+bc3_length+sp2_length+bc2_length
	bc1_start = UMI_length+bc3_length+sp2_length+bc2_length+sp1_length
	bc1_end = UMI_length+bc3_length+sp2_length+bc2_length+sp1_length+bc1_length
		
	line1 = f1.readline()
	line2 = f2.readline()
	count = 0
	while (line1):
		count +=1
		if count % 10000 == 0:
			print ("We have processed %s number of reads"%(count))
		name_r1=line1
		name_r2=line2
		line1 = f1.readline()
		line2 = f2.readline()
		bc3 = line1[bc3_start:bc3_end]
		bc2 = line1[bc2_start:bc2_end]
		bc1 = line1[bc1_start:bc1_end]
		UMI = line1[:UMI_length]
		bc1_dist = find_dist(bc1_kmer,bc1,barcode_1_list)
		bc2_dist = find_dist(bc2_kmer,bc2,barcode_2_list)
		bc3_dist = find_dist(bc3_kmer,bc3,barcode_3_list)

		if bc1_dist and bc2_dist and bc3_dist :
			barcode_dict[bc3][bc2][bc1]+=1
			first_line_r1 = '@' + ",".join([UMI,bc3,bc2,bc1]) + ',' + name_r1[1:]
			first_line_r2 = '@' + ",".join([UMI,bc3,bc2,bc1])+ ',' + name_r2[1:]
			matched_list_R1.append(first_line_r1)
			matched_list_R2.append(first_line_r2)

			matched_list_R1.append(line1)
			matched_list_R2.append(line2)	
			
			third_line_r1 = f1.readline()
			third_line_r2 = f2.readline()
			matched_list_R1.append(third_line_r1)
			matched_list_R2.append(third_line_r2)	
	
			four_line_r1 = f1.readline()
			four_line_r2 = f2.readline()
			matched_list_R1.append(four_line_r1)
			matched_list_R2.append(four_line_r2)	
	

		else:
			junk_list_R1.append(name_r1)
			junk_list_R2.append(name_r2)
			junk_list_R1.append(line1)
			junk_list_R2.append(line2)
			
			third_line_r1 = f1.readline()
			third_line_r2 = f2.readline()
			junk_list_R1.append(third_line_r1)
			junk_list_R2.append(third_line_r2)	
	
			four_line_r1 = f1.readline()
			four_line_r2 = f2.readline()
			junk_list_R1.append(four_line_r1)
			junk_list_R2.append(four_line_r2)	
				
			
		line1 = f1.readline()
		line2 = f2.readline()

	junk_R1_output = output_folder + "/" + label + ".junk.R1.fastq.gz"
	junk_R2_output = output_folder + "/" + label + ".junk.R2.fastq.gz"
	output_to_fastq_gz(junk_R1_output,junk_list_R1)
	output_to_fastq_gz(junk_R2_output,junk_list_R2)
	
	matched_R1_output = output_folder + "/" + label + ".matched.R1.fastq.gz"
	matched_R2_output = output_folder + "/" + label + ".matched.R2.fastq.gz"
	output_to_fastq_gz(matched_R1_output,matched_list_R1)
	output_to_fastq_gz(matched_R2_output,matched_list_R2)
	df = dict3d_to_df(barcode_dict)
	df.to_csv(output_folder + "/" + label + ".total_number_reads.tsv",sep="\t",index=False)
	# for bc3 in barcode_3_list:
		# for bc2 in barcode_2_list:
			# for bc1 in barcode_1_list:
				# lines = matched_list_R1
				# if len(lines) < 4:
					# continue
				# output_R1 = output_folder + "/" + label + ".%s.%s.%s.R1.fastq.gz"%(bc3,bc2,bc1)
				# output_R2 = output_folder + "/" + label + ".%s.%s.%s.R2.fastq.gz"%(bc3,bc2,bc1)
				# output_to_fastq_gz(output_R1,matched_list_R1)
				# output_to_fastq_gz(output_R2,matched_list_R2)
			
def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	mainParser.add_argument("-r1","--read1",  help="R1 fastq file", required=True)
	mainParser.add_argument("-r2","--read2",  help="R2 fastq file", required=True)
	mainParser.add_argument("--sample_ID",  help="sample_ID", required=True)
	# mainParser.add_argument("--RT_primer",  help="RT_primer", required=True)
	mainParser.add_argument("--barcode_1_list",  help="list of barcodes 1", required=True)
	mainParser.add_argument("--barcode_2_list",  help="list of barcodes 2", required=True)
	mainParser.add_argument("--barcode_3_list",  help="list of barcodes 3", required=True)
	mainParser.add_argument("--UMI_length",  help="barcode 3 allowed number of mismatches", type=int, required=True)
	mainParser.add_argument("--bc3_length",  help="barcode 3 allowed number of mismatches", type=int, required=True)
	mainParser.add_argument("--sp2_length",  help="barcode 3 allowed number of mismatches", type=int, required=True)
	mainParser.add_argument("--bc2_length",  help="barcode 3 allowed number of mismatches", type=int, required=True)
	mainParser.add_argument("--sp1_length",  help="barcode 3 allowed number of mismatches", type=int, required=True)
	mainParser.add_argument("--bc1_length",  help="barcode 3 allowed number of mismatches", type=int, required=True)
	mainParser.add_argument("--BC1_error",  help="barcode 1 allowed number of mismatches", type=int,default=0)
	mainParser.add_argument("--BC2_error",  help="barcode 2 allowed number of mismatches", type=int,default=0)
	mainParser.add_argument("--BC3_error",  help="barcode 3 allowed number of mismatches", type=int,default=0)
	

	mainParser.add_argument('--input', help=argparse.SUPPRESS)

	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def main():

	args = my_args()
	barcode_1_list = pd.read_csv(args.barcode_1_list,header=None)[0].tolist()
	barcode_2_list = pd.read_csv(args.barcode_2_list,header=None)[0].tolist()
	barcode_3_list = pd.read_csv(args.barcode_3_list,header=None)[0].tolist()

	
	sci_l3_demultiplexing(args.read1,args.read2,args.sample_ID, barcode_1_list, barcode_2_list, barcode_3_list, args.BC1_error, args.BC2_error, args.BC3_error,args.UMI_length,args.bc3_length,args.sp2_length,args.bc2_length,args.sp1_length,args.bc1_length)
	

if __name__ == "__main__":
	main()
	
