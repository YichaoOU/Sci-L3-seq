# Sci-L3-seq

Demultiplexing, read mapping, QC for sci-L3-seq, WGS, target-seq, and co-assay

Reference code was modified and fit into HemTools pipeline (https://hemtools.readthedocs.io/en/latest/content/Gallery/run_lsf.html).

# Input

Copy (or use `ln -s`) your fastq files to your working directory, prepare fastq.tsv using `run_lsf.py --guess_input`.

The input file `fastq.tsv` looks like below:

	Banana.R1.fastq.gz	Banana.R2.fastq.gz	Banana
	Orange.R1.fastq.gz	Orange.R2.fastq.gz	Orange

# Usage for example data

The following is supposed to work at St. Jude HPC.

```

hpcf_interactive

export PATH=$PATH:"/home/yli11/HemTools/bin"

git clone https://github.com/YichaoOU/Sci-L3-seq

cd Sci-L3-seq

vim pipeline/sci_l3_seq.lsf

# find and replace src=/home/yli11/Programs/Sci-L3-seq/src with your correct abs path

cd data

module load python/2.7.13

run_lsf.py --guess_input

run_lsf.py -f fastq.tsv --user_lsf ../pipeline/sci_l3_seq.lsf -j my_expamle_results

```

If you can't figure out how to run it in your local HPC, please feel free to contact me.

# Output

sample QC (e.g., number of reads contains all the barcodes, collision rate) summary is provided in `jid/sample_QC.tsv`

An example of collision plot:

![example](collision_plot.png)

# Notes

1. P5 or P7 sequencing primer are equally possible to be added at the UMI end or gRNA end. So de-multiplexing should look at both R1 and R2 reads. 

2. The reads that we can used to do demultiplexing (R1 or R2) should have this format: UMI (4nt) + SSS_barcode (6nt) + GGGATGCAGCTCGCTCCTG (20nt, RT_primer) + barcode_2 (7nt) + spacer_sequence (6nt, GTCTTG) + barcode_1 (8nt) + Tn5 (19nt, AGATGTGTATAAGAGACAG) + gRNA

3. By default, barcode 3 allows no mismatch, barcode 2 and barcode 1 each allows 1 mismatch, RT primer allows 3 mismatch.

# fastq read example

![example](fastq_read_example.png)

The read structure helps to follow the pipeline described below.



# Pipeline

step1:

	1. given RT_primer, assign PE reads to junk (noRT.fastq.gz) and not_junk

	2. swap R2 R1 if R1 doesn't have RT_primer

	3. match barcode_3, not matched reads will be discarded (noBC3.fastq.gz)

	4. output barcode_3-RT PE reads to R1.ordered.fastq.gz and R2.ordered.fastq.gz

step2:

	5. given barcode_1, barcode_2, barcode_3, parse files from step 4 to matched or junk, renamed read name if matched
	
step3:

	6. cutadapt trim using Tn5 sequence 
	
step4:

	7. bwa mapping, 
	
step5:

	8. bedtools bamtobed
	
step 6:

	9. summerize results to table and figure, provided step3_QC_summary.py  step4_calculate_collision_rate.py




## Reference

https://github.com/Yue-Jiang/sciliantifig


