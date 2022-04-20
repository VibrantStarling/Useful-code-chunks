#!/bin/bash

# This chunk is only suited to single paired end read samples that already have their own directory
# Please adapt this code to your own filing system if you have multiple samples to process!

# Read in from an input text file of name prefixes
for i in $(cat $1)
do

# Run phyloFlash OPTION1 (sequential)
# Get first 10000 read from fastq file and compute average read length
N=`zcat ${i}_1.fastq.gz | head -n 10000 | awk 'NR%4 == 2 { count++; bases += length($0)} END {printf "%3.0f\n", bases/count}'`
echo 'average read length in' ${i}': ' ${N}

#run phyloflash (change the pathway to phyloflash.pl if it is saved dfferently on your system
~/phyloFlash-pf3.4/phyloFlash.pl -lib ${i} -dbhome ~/132 -read1 ${i}_1.fastq.gz -read2 ${i}_2.fastq.gz -zip -CPUs 48 -readlength ${N} -taxlevel 5

done
