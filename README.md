# Useful-code-chunks
Pieces of code that make life easier, but I would otherwise probably forget the day after I have worked it out

# File handling
1. [Mass rename files in different directories](#1)
2. [Replace multiple pairs of patterns from a tab delimited list](#10)

# Running programmes with additional bits
1. [Run PhyloFlash and get average read lengths](#2)
2. [Process and align paired rnaseq reads](#12)

# Genome statistics
1. [Get No. bp for multiple files](#3)
2. [Get read coverage for a list of contigs](#4)
3. [Calculate GC content](#5)
4. [Calculate N50, N90, length, number of contigs and largest contig](#6)
5. [Calculate average CDS size](#7)
9. [Get aa sequence from a gff](#8)
10. [Calculate reciprocal best hits](#9)
11. [Get bam file statistics with samtools](#11)

## Mass rename files in different directories <a name="1"></a>
> How to mass rename files with specific prefixes in multiple different directories so file names are more informative. (Just incase you forget to add the option in programmes like metabat and have about 500 files that need prefixes)

```
# ADDING PREFIXES METABAT BINS IN DIFFERENT FOLDERS
# needs a list of directory names/prefixes you want to add to files further in. 
# E.g directory for SRA data named DRR215885 with metabat bins inside called bin.1.fa bin.2.fa ... etc. 
# This adds DRR215885 as a prefix to each bin to become DRR215885_bin.1.fa, DRR215885_bin.2.fa ... etc.

for i in $(cat names); do cd ${i}/*metabat*/; for file in $(ls *.fa); do  name="${i}_${file}"; mv -vi ./${file} ./${name};  done; cd ../../; done
```

## Replace multiple pairs of patterns from a tab delimited list <a name="10"></a>
> Do a whole bunch of replacements on one file from a tab delimited list
```
while read a b; do sed -i "s/${b}/${a}/g" target-file.txt; echo ${a}; done < pattern-pairs-list.txt
```
The `pattern-pairs-list.txt` should be formatted like:
```
new_pattern_1  old_pattern_1
new_pattern_2  old_pattern_2
new_pattern_3  old_pattern_3
new_pattern_4  old_pattern_4
```

## Running Phyloflash and computing average read length <a name="2"></a>
> Calculate your read length and feed it directly into [PhyloFlash](http://hrgv.github.io/phyloFlash/ "PhyloFlash Manual"). 
> This chunk is only suited to single paired end read samples that already have their own directory. Please adapt this code to your own filing system if you have multiple samples to process or if you save your scripts in a different way! Also check your fastq file names match the format in this code chunk!
```
# run runPhloflash.sh with a text list of file prefixes

bash ~/scripts/runPhyloflash.sh [inputfile name]
```
## Process and align paired rnaseq reads <a name="12"></a>
This script will softmask a genome to prepare it for rnaseq alignment. It will then take a list of SRA numbers for paired rnaseq data, trim those files and align them to the masked genome.

**Dependancies**
```
# singularity image for dfam-TEtool (dfam-tetools-latest.sif) in your home directory
singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
# hisat2
conda install bioconda::hisat2
# trimmomatic
conda install bioconda::trimmomatic
# samtools
conda install bioconda::samtools
```
**Inputs**

1 - a list of SRA numbers (e.g. SRR6493555) for paired rnaseq sets with each number on a new line

2 - the genome to mask and align rnaseq to. Must end in the extension `.fna`

3 - the name for the masked genome database

**Useage**
```
bash process_and_align_paired_rnaseq.sh sra_list genome.fna genome_db_name_prefix
```

## get the number of bp for a bunch of fasta files <a name="3"></a>

```
for i in $(ls *.fa)
do 
grep -v ">" ${i} | wc | awk -v var=${i} '{print var," = ",$3-$1}'
done
```

## Get coverage for a metabat bin from megahit final contigs depths <a name="4"></a>
```
## get a list of contig names
grep ">" bin.fa | sed 's/>//g' > contigs.txt

## cross reference the depth file from the megahit assembly
for i in $(cat contigs.txt)
do grep -w ${i} contigs.fa.depth.txt >> bin-depths.txt
done

## get average bin depth
cut -f3 bin-depths.txt | awk '{s+=$1}END{print "ave:",s/NR}'
```
## Calculate GC content <a name="5"></a>
use `gc-calculator.sh` as follows
```
# for one file
bash ~/scripts/gc-calculator.sh fasta.fa

# for multiple files
for i in $(ls *.fa); do bash ~/scripts/gc-calculator.sh ${i}; done
```
## Calculate N50, N90, length, number of contigs and largest contig <a name="6"></a>
Use this script https://github.com/hcdenbakker/N50.sh/blob/682bc724ecad22559e2fcd8ef8bcfc48ed1e8e5f/N50.sh
```
bash N50-calculator.sh your-genome.fa
```

## Calculate average CDS size <a name="7"></a>
you will need a prokka ouput .tsv file as an input
```
cat *.tsv | awk '$2 == "CDS"'| awk '{total += $3} END {print total/NR}'
```

## Calculate coding density <a name="8"></a>
Prokka will give you the number of predicted CDS, and you can use its tsv output to calculate the average CDS length.
`$1 = avg. CDS length (bp)`
`$2 = number of CDS`
`$3 = genome length (bp)`
```
awk 'BEGIN {print ($1*$2/$3*100)}'
```

## Get aa sequence for the an exon (mRNA) from a gff and fasta <a name="9"></a>
AGAT has a whole bunch of very useful tools for getting information out of and editing gff files, including this one:
https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html

This is how to get the mRNA aa sequence:
```
agat_sp_extract_sequences.pl -g GFF -f corresponding_fasta -t exon --aa --merge -o output_name.faa
```

## Calculate reciprocal best hits with `find-reciprocal-best-hits.py` <a name="10"></a>

Calculate reciprocal best hits and make some nice histograms and density plots
Based on code from [here](https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html)

**Usage:**
```
find-reciprocal-best-hits.py [-h] -s1 QUERY_AA -s2 SUBJECT_AA -o OUTPUT

optional arguments:
  -h, --help    show this help message and exit
  -s1 QUERY_AA, --query_aa QUERY_AA   
                one set of aa sequence
  -s2 SUBJECT_AA, --subject_aa SUBJECT_AA    
                second set of aa sequence
  -o OUTPUT, --output OUTPUT    
                Name for outputs
```
**Output:**
Blast results for forward and reverse reads:

-`outputName-fwd-results.tab`
-`outputName-rev-results.tab`

Reciprocal best hits:

-`Cposadasi-test-rbbh.csv`
- histograms of normalised bitscores for forward, reverse and reciprocal reads:
- 2D density plots of normalised bitscores for forward, reverse and reciprocal reads:

## Get bam file statistics with samtools <a name="11"></a>
Run this in a directory that contains bam files, preferably sorted bam files

```
for i in $(ls *.bam)
do echo "SAMTOOLS FLAGSTAT" > ${i}_stats
echo "" >> ${i}_stats
samtools flagstat ${i} >> ${i}_stats
echo "" >> ${i}_stats
echo "" >> ${i}_stats
echo "SAMTOOLS STAT" >> ${i}_stats
echo "" >> ${i}_stats
samtools stats ${i} |grep ^SN | cut -f 2- >> ${i}_stats
done
```
