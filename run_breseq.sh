#!/bin/bash
#SBATCH --time=10:00:00 
#SBATCH --partition=gelifes 
#SBATCH --mail-type=ALL
#SBATCH --job-name=breseq_ABCD_D7_w_filters
#SBATCH --output=breseq_ABCD_D7_w_filters-%j.out
#SBATCH --mem=128000

#load module
module load breseq/0.34.1a-foss-2018a

#create array where each row is a different sample read (either forward or reverse)
readarray -t R1 < <(find -name "*R1*.gz")
readarray -t R2 < <(find -name "*R2*.gz")

#set IFS and sort
IFS=$'\n' 
sorted_R1=($(sort -d <<<"${R1[*]}"))
sorted_R2=($(sort -d <<<"${R2[*]}"))

unset IFS

#pipe loop to parallel command
#set breseq parameters
cnt=${#sorted_R2[@]} 
for ((i=0;i<cnt;i++)); do
breseq -p \
--polymorphism-minimum-variant-coverage 4 \
--polymorphism-minimum-total-coverage 10 \
--polymorphism-bias-cutoff 0.05 \
-r A-Ancestral_edited.gbk \
-r B-Ancestral_edited.gbk \
-r C-Ancestral_edited.gbk \
-r D-Ancestral_edited.gbk \
-o "${sorted_R1[i]%%_3*}" \
${sorted_R1[i]} \
${sorted_R2[i]} \
&> log.txt

done \
| parallel -P 16

#make function out of breseq command
func breseq_run()
{	breseq -p \
	--polymorphism-minimum-variant-coverage 4 \
	--polymorphism-minimum-total-coverage 10 \
	--polymorphism-bias-cutoff 0.05 \
	-r A-Ancestral_edited.gbk \
	-r B-Ancestral_edited.gbk \
	-r C-Ancestral_edited.gbk \
	-r D-Ancestral_edited.gbk \
	-o ${1%%_3*} \
	$1 \
	$2 \
	&> log.txt
}
export -f breseq_run

#run function in parallel
parallel --jobs 16 breseq_run




