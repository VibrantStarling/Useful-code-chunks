#! /bin/bash

#Get the file name

if [ $# != 1 ]; then
        echo "Please include a fasta input on the command line!"
        exit
fi


#Find all Gs and Cs, sequence length and then calculate GC content
lengths=`grep -v ">" $1 | wc | awk '{print $3-$1}'`
C=`grep -v ">" $1 | awk -F"C" '{print NF-1}'| awk '{sum += $1} END {print sum}'`
G=`grep -v ">" $1| awk -F"G" '{print NF-1}'| awk '{sum += $1} END {print sum}'`
GC=`expr $C + $G`
percent=`echo "$GC / $lengths * 100" | bc -l`

#Print result

echo "GC Content for $1 is  $percent"
