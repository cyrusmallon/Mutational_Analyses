# Variant Discovery in Bacterial Genomes
This reposity details the processing of high throughput genomic data, from sequences to (eventually) visualization. 

## run_breseq.sh
This is a shell script for  passing the quality controlled sequencing data (i.e., the short read data coming from the sequencer, in this case illumina data) through the breseq pipeline. This script is designed to run n number of read samples in parallel to a computer cluster running SLURM scheduling. Because I am using mixed communities of bacteria from one to four species, there are four reference genomes in this particualar script upon which the short reads are mapped. 
