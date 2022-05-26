# Useful-code-chunks
Pieces of code that make life easier, but I would otherwise probably forget the day after I have worked it out

## File juggling
> How to mass rename files with specific prefixes in multiple different directories so file names are more informative. (Just incase you forget to add the option in programmes like metabat and have about 500 files that need prefixes)

```
# ADDING PREFIXES METABAT BINS IN DIFFERENT FOLDERS
# needs a list of directory names/prefixes you want to add to files further in. 
# E.g directory for SRA data named DRR215885 with metabat bins inside called bin.1.fa bin.2.fa ... etc. 
# This adds DRR215885 as a prefix to each bin to become DRR215885_bin.1.fa, DRR215885_bin.2.fa ... etc.

for i in $(cat names); do cd ${i}/*metabat*/; for file in $(ls *.fa); do  name="${i}_${file}"; mv -vi ./${file} ./${name};  done; cd ../../; done
```


## Running Phyloflash and computing average read length
> Calculate your read length and feed it directly into [PhyloFlash](http://hrgv.github.io/phyloFlash/ "PhyloFlash Manual"). 
> This chunk is only suited to single paired end read samples that already have their own directory. Please adapt this code to your own filing system if you have multiple samples to process or if you save your scripts in a different way! Also check your fastq file names match the format in this code chunk!
```
# run runPhloflash.sh with a text list of file prefixes

bash ~/scripts/runPhyloflash.sh [inputfile name]
```
## get the number of bp for a bunch of fasta files

```
for i in $(ls *.fa)
do 
grep -v ">" ${i} | wc | awk -v var=${i} '{print var," = ",$3-$1}'
done
```

## Get coverage for a metabat bin from megahit final contigs depths
```
## get a list of contig names
grep ">" bin.fa | sed 's/>//g' > contigs.txt

## cross reference the deopth file from the megahit assembly
for i in $(cat contigs.txt)
do grep -w ${i} contigs.fa.depth.txt >> bin-depths.txt
done

## get average bin depth
cut -f3 bin-depths.txt | awk '{s+=$1}END{print "ave:",s/NR}'
```
## Calculate GC content
use `gc-calculator.sh` as follows
```
# for one file
bash ~/scripts/gc-calculator.sh fasta.fa

# for multiple files
for i in $(ls *.fa); do bash ~/scripts/gc-calculator.sh ${i}; done
```
## Calculate N50, N90, length, number of contigs and largest contig
Use this script https://github.com/hcdenbakker/N50.sh/blob/682bc724ecad22559e2fcd8ef8bcfc48ed1e8e5f/N50.sh
