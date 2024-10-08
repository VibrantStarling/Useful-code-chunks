#!/bin/bash

# adapted from https://www.biostars.org/p/139006/#139030
# takes multiple fastq files and calculates the number of reads in them

if [ $# -lt 1 ] ; then
    echo ""
    echo "usage: count_fastq.sh [fastq_file1] <fastq_file2> ..|| *.fastq"
    echo "counts the number of reads in a fastq file"
    echo ""
    exit 0
fi

filear=${@};
for i in ${filear[@]}
do
counts=$(pigz -dc ${i} | awk 'NR % 4 == 2' | wc -l)
echo -n -e "\t$i : "
echo $counts | \
sed -r '
  :L
  s=([0-9]+)([0-9]{3})=\1,\2=
  t L'
done
