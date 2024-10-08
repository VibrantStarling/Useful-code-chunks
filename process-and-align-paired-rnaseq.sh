#! /bin/bash
# This script will softmask a genome to prepare it for rnaseq alignment. 
# It will then take a list of SRA numbers for paired rnaseq data, trim those 
# files and align them to the masked genome.

# DEPENDANCIES
# singularity image for dfam-TEtool (dfam-tetools-latest.sif) in your home directory
#        Use: singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
# hisat2
#        Use: conda install bioconda::hisat2
# trimmomatic
#        Use: conda install bioconda::trimmomatic
# samtools
#        Use: conda install bioconda::samtools
# fastq-pair
#       conda install bioconda::fastq-pair

# INPUTS
# 1 - a list of SRA numbers (e.g. SRR6493555) for paired rnaseq sets with each number on a new line
# 2 - the genome to mask and align rnaseq to
# 3 - the name for the masked genome database

 if [ ! $# -eq 3 ]
    then
    echo " Three arguments needed:  "
    echo " the sra list"
    echo " GENOME.fna to align rnaseq to" 
    echo " A name to give the genome database" 
    exit 1
 fi

SRA_LIST=$1
GENOME=$2
DB=$3
IDX=$(basename -s .fna ${GENOME})

# soft mask reads with TETools repeatmodler2 and repeatmasker. 
#    > can struggle if memory is not fast enough and data may need to be stored on compute HPC nodes
#    > This may need to be done more rigourously for larger genomes
tput setaf 5; echo "softmasking genome"; tput sgr0
~/dfam-tetools.sh
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/dfam-tetools-latest.sif BuildDatabase -name ${DB} ${GENOME}
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/dfam-tetools-latest.sif RepeatModeler -database ${DB} -threads 32 -LTRStruct
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/dfam-tetools-latest.sif RepeatMasker -pa 32 -lib ${DB}-families.fa -xsmall ${GENOME}
tput setaf 2; echo "softmasking complete"; tput sgr0


# clean and align SRA to the genome
for RNA_PREFIX in $(cat ${SRA_LIST})
do
    # define the file objects
    RNASEQ_FWD=${RNA_PREFIX}_1.fastq.gz
    RNASEQ_REV=${RNA_PREFIX}_2.fastq.gz
    
    # trim the fastq file with trimmomatic
    tput setaf 6; echo "------START of trimming for ${RNA_PREFIX}"; tput sgr0
    trimmomatic PE -phred33 -threads 32 ${RNASEQ_FWD} ${RNASEQ_REV} ${RNA_PREFIX}_fpaired.fq.gz ${RNA_PREFIX}_funpaired.fq.gz ${RNA_PREFIX}_rpaired.fq.gz ${RNA_PREFIX}_runpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
    echo "Number of trimmed forward paired reads: " 
    echo $(zcat ${RNA_PREFIX}_rpaired.fq.gz |wc -l)/4|bc
    echo "Number of trimmed reverse paired reads: " 
    echo $(zcat ${RNA_PREFIX}_fpaired.fq.gz |wc -l)/4|bc
    tput setaf 2; echo "------END of  trimming for ${RNA_PREFIX}------"; tput sgr0
    echo

    # define the trimmed file objects
    RNASEQ_FWD=${RNA_PREFIX}_fpaired.fq.gz
    RNASEQ_REV=${RNA_PREFIX}_rpaired.fq.gz

    # paired with a list of SRA names (SRRXXXXXX) to align 
    tput setaf 6; echo "------START of alignment for ${RNA_PREFIX}"; tput sgr0
    time hisat2 -p 32 -q -x ${IDX} -1 ${RNASEQ_FWD} -2 ${RNASEQ_REV} > ${RNA_PREFIX}-hisat2-rnaseq.sam  2> ${RNA_PREFIX}-hisat2-align.err
    samtools view -@ 12 -bS ${RNA_PREFIX}-hisat2-rnaseq.sam -o ${RNA_PREFIX}-hisat2-rnaseq.bam
    samtools sort -@ 12 ${RNA_PREFIX}-hisat2-rnaseq.bam -o ${RNA_PREFIX}-hisat2-rnaseq_sorted.bam
    rm ${RNA_PREFIX}-hisat2-rnaseq.sam; rm ${RNA_PREFIX}-hisat2-rnaseq.bam
    tput setaf 2; echo "------END of  alignment for ${RNA_PREFIX}------"; tput sgr0

done
