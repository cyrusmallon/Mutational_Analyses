# Variant Discovery in Bacterial Genomes
This reposity details the processing of high throughput genomic data, from sequences to visualization. 

## run_breseq.sh
This is a shell script for  passing the quality controlled sequencing data (i.e., the short read data coming from the sequencer, in this case illumina data) through the breseq pipeline. This script is designed to run n number of read samples in parallel to a computer cluster running SLURM scheduling. Because I am using mixed communities of bacteria from one to four species, there are four reference genomes in this particualar script upon which the short reads are mapped. 

## summarise_breseq_results.py
This is a python wrapper/parser script that summarises all breseq results into a single .csv file. Simply create a directory with all breseq output folders and this script inside the parent directy and run the script to create a .csv file. This script was originally written by Ali May (linkedin.com/in/alimay) and later debugged (minor changes) by myself and Thijs Janzen (thijsjanzen@gmail.com). There is still an issue with mutational calls that span two rows, as they are called a single mutation in this script. Feel free to fix.
