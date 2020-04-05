#!/usr/bin/env python

import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import os
import argparse

'''
	this script accepts and opens a R1 file, a R2 file, an RT primer sequence, and a SSS barcode list; then reorders the read pairs, and extracts and matches the nearest sss barcode and attach the sss barcode to R1 and R2 read names.
	the output files include ordered R1 and R2 files, noRTprimer.R1.fq.gz and noRTprimer.R2.fq.gz
	
'''

def SSS_barcode_attach_list(Read1,Read2,label, barcode_list, RT_primer, mismatch_sss, mismatch_RT):
	#open the read1, read2, and output files
	output_folder = label
	R1_out_file = output_folder + "/" + label + ".R1.ordered.fastq.gz"
	R2_out_file = output_folder + "/" + label + ".R2.ordered.fastq.gz"
	noRT_R1_out_file = output_folder + "/" + label + ".noRTprimer.R1.fq.gz"
	noRT_R2_out_file = output_folder + "/" + label + ".noRTprimer.R2.fq.gz"

	f1 = gzip.open(Read1)
	f2 = gzip.open(Read2)
	f3 = gzip.open(R1_out_file, 'wb')
	f4 = gzip.open(R2_out_file, 'wb')
	f5 = gzip.open(noRT_R1_out_file, 'wb')
	f6 = gzip.open(noRT_R2_out_file, 'wb')
	
	line1 = f1.readline()
	line2 = f2.readline()

	while (line1):
		name_r1=line1
		name_r2=line2
		line1 = f1.readline()
		line2 = f2.readline()
		target_rt_r1 = line1[9:29]
		target_rt_r2 = line2[9:29]
		target_sss_r1 = line1[4:10]
		target_sss_r2 = line2[4:10]
		
		junk=True
		dist_rt_r1 = distance(RT_primer, target_rt_r1)
		dist_rt_r2 = distance(RT_primer, target_rt_r2)
		
		if (dist_rt_r1 <= mismatch_RT):
			junk=False
			for barcode in barcode_list:
				mismatch = distance(barcode, target_sss_r1)
				find = False
				
				if (mismatch <= mismatch_sss):
					find = True
					UMI = line1[:4]
					first_line_r1 = '@' + barcode + ',' + UMI + ',' + name_r1[1:]
					first_line_r2 = '@' + barcode + ',' + UMI + ',' + name_r2[1:]
					f3.write(first_line_r1)
					f4.write(first_line_r2)
	
					second_line_r1 = line1
					second_line_r2 = line2
					f3.write(second_line_r1)
					f4.write(second_line_r2)
	
					third_line_r1 = f1.readline()
					third_line_r2 = f2.readline()
					f3.write(third_line_r1)
					f4.write(third_line_r2)
	
					four_line_r1 = f1.readline()
					four_line_r2 = f2.readline()
					f3.write(four_line_r1)
					f4.write(four_line_r2)
	
					line1 = f1.readline()
					line2 = f2.readline()
					break
					
			if find == False:
				line1 = f1.readline()
				line1 = f1.readline()
				line1 = f1.readline()
				line2 = f2.readline()
				line2 = f2.readline()
				line2 = f2.readline()
				
		#swap R1 R2 but keep read name
		elif (dist_rt_r2 <= mismatch_RT):
			junk=False
			for barcode in barcode_list:
				mismatch = distance(barcode, target_sss_r2)
				find = False
				
				if (mismatch <= mismatch_sss):
					find = True
					UMI = line2[:4]
					first_line_r1 = '@' + barcode + ',' + UMI + ',' + 'swap' + ',' + name_r1[1:]
					first_line_r2 = '@' + barcode + ',' + UMI + ',' + 'swap' + ',' + name_r2[1:]
					f3.write(first_line_r1)
					f4.write(first_line_r2)
	
					second_line_r1 = line2
					second_line_r2 = line1
					f3.write(second_line_r1)
					f4.write(second_line_r2)
	
					third_line_r1 = f2.readline()
					third_line_r2 = f1.readline()
					f3.write(third_line_r1)
					f4.write(third_line_r2)
	
					four_line_r1 = f2.readline()
					four_line_r2 = f1.readline()
					f3.write(four_line_r1)
					f4.write(four_line_r2)
	
					line1 = f1.readline()
					line2 = f2.readline()
					break
					
			if find == False:
				line1 = f1.readline()
				line1 = f1.readline()
				line1 = f1.readline()
				line2 = f2.readline()
				line2 = f2.readline()
				line2 = f2.readline()
				
		#write junk file
		else: 
			junk = True
			first_line_r1 = name_r1
			first_line_r2 = name_r2
			f5.write(first_line_r1)
			f6.write(first_line_r2)

			second_line_r1 = line1
			second_line_r2 = line2
			f5.write(second_line_r1)
			f6.write(second_line_r2)

			third_line_r1 = f1.readline()
			third_line_r2 = f2.readline()
			f5.write(third_line_r1)
			f6.write(third_line_r2)

			four_line_r1 = f1.readline()
			four_line_r2 = f2.readline()
			f5.write(four_line_r1)
			f6.write(four_line_r2)
			line1 = f1.readline()
			line2 = f2.readline()
			
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()

def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# mainParser.add_argument('-j',"--jid",  help="enter a job ID, which is used to make a new directory. Every output will be moved into this folder.", default=current_file_base_name+'_'+username+"_"+str(datetime.date.today()))	
	# mainParser.add_argument("--bamCoverage_addon",  help="for PE data, you add --center to get sharper peaks", default="")
	mainParser.add_argument("-r1","--read1",  help="R1 fastq file", required=True)
	mainParser.add_argument("-r2","--read2",  help="R2 fastq file", required=True)
	mainParser.add_argument("--sample_ID",  help="sample_ID", required=True)
	mainParser.add_argument("--RT_primer",  help="RT_primer", required=True)
	mainParser.add_argument("--barcode_list",  help="list of barcodes", required=True)
	mainParser.add_argument("--mismatch_sss",  help="mismatch_sss", type=int,default=0)
	mainParser.add_argument("--mismatch_RT",  help="mismatch_RT", type=int,default=3)
	# mainParser.add_argument("--barcode_list",  help="list of barcodes", required=True)
	
	mainParser.add_argument('--input', help=argparse.SUPPRESS)

	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def main():

	args = my_args()
	SSS_barcode_attach_list(args.read1, args.read2, args.sample_ID, args.barcode_list, args.RT_primer, args.mismatch_sss, args.mismatch_RT)
	

if __name__ == "__main__":
	main()
	
