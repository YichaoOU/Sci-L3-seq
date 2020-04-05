# Sci-L3-seq
Demultiplexing code for sci-L3-seq, WGS, target-seq, and co-assay

Reference code was modified and fit into HemTools pipeline (https://hemtools.readthedocs.io/en/latest/content/Gallery/run_lsf.html).



# Methods

1. P5 or P7 sequencing primer are equally possible to be added at the UMI end or gRNA end. So de-multiplexing should look at both R1 and R2 reads. 

2. The reads that we can used to do demultiplexing (R1 or R2) should have this format: UMI (4nt) + SSS_barcode (6nt) + GGGATGCAGCTCGCTCCTG (20nt, RT_primer) +


## Reference

https://github.com/Yue-Jiang/sciliantifig


