#!/bin/bash
#Script to perform a reciprocal blast search
#From https://morphoscape.wordpress.com/2020/08/18/reciprocal-best-hits-blast-rbhb/
#Usage: bash runRBLAST.sh PATH/TO/QUERY/FILE PATH/TO/DB/FILE PATH/TO/OUTPUTS
#Usage ex: bash runRBLAST.sh PATH/TO/INPUT1/species1.fasta.transdecoder.pep PATH/TO/INPUT2/species2.fasta.transdecoder.pep PATH/TO/OUTPUTS

#Input query file
inputQuery="$1"
#Input DB reciprocal file
inputDB="$2"
#Path to output results
outputPath="$3"


#Move to query directory
queryPath=$(dirname "$inputQuery")
cd "$queryPath"
#Make blastable protein DB of the input query
makeblastdb -in "$inputQuery" -dbtype prot

#Move to DB directory
dbPath=$(dirname "$inputDB")
cd "$dbPath"
#Make blastable protein DB of the input DB
makeblastdb -in "$inputDB" -dbtype prot

#Output start status message
echo "Beginning reciprocal BLAST..."

#Move to outputs folder
mkdir "$outputPath"
cd "$outputPath"

#Use blastp to search a database
blastp -query "$inputQuery" -db "$inputDB" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > fwd_blast.tsv

#Switch query and search paths for reciprocal search
blastp -query "$inputDB" -db "$inputQuery" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > reverse_blast.tsv

#Output end status message
echo "Finished reciprocal BLAST!"