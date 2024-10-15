# BRAKER3 pipeline

'''
# have seqkit installed for fastq file processing to keep trimmomatic happy

# How to run all versions of BRAKER and GALBA (feat. OMArk and BUSCO)
https://www.youtube.com/watch?v=UXTkJ4mUkyg

# How to access BRAKER3:
https://hub.docker.com/r/teambraker/braker3

singularity build braker3.sif docker://teambraker/braker3:latest
singularity exec braker3.sif braker.pl
export BRAKER_SIF=braker3.sif
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
bash test1.sh # tests BRAKER1
bash test2.sh # tests BRAKER2
bash test3.sh # tests BRAKER3


# How to access TETools:
https://github.com/Dfam-consortium/TETools
# make sure perl is installed through perlbrew, NOT conda
# make sure h5py is installed with pip

curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh
./dfam-tetools.sh

'''
GENOME=""
DB=""
IDX=""

# trim paired end reads with trimmomatic
# replace TruSeq3-PE.fa with the path to whatever adapter file you're using
SRA_LIST=sra-list.txt
for RNA_PREFIX in $(cat ${SRA_LIST})
do
    # define the file objects
    RNASEQ_FWD=${RNA_PREFIX}_1.fastq.gz
    RNASEQ_REV=${RNA_PREFIX}_2.fastq.gz
    
    # remove spaces from the fastq
    seqkit replace -p "\s.+" ${RNASEQ_FWD} -o ${RNA_PREFIX}_1_clean.fastq.gz
    seqkit replace -p "\s.+" ${RNASEQ_REV} -o ${RNA_PREFIX}_2_clean.fastq.gz

    RNASEQ_FWD=${RNA_PREFIX}_1_clean.fastq.gz
    RNASEQ_REV=${RNA_PREFIX}_2_clean.fastq.gz
    
    # trim the fastq file with trimmomatic
    tput setaf 6; echo "------START of trimming for ${RNA_PREFIX}"; tput sgr0
    trimmomatic PE -phred33 -threads 32 ${RNASEQ_FWD} ${RNASEQ_REV} ${RNA_PREFIX}_fpaired.fq.gz ${RNA_PREFIX}_funpaired.fq.gz ${RNA_PREFIX}_rpaired.fq.gz ${RNA_PREFIX}_runpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
    echo "Number of trimmed forward paired reads: " 
    echo $(zcat ${RNA_PREFIX}_rpaired.fq.gz |wc -l)/4|bc
    echo "Number of trimmed reverse paired reads: " 
    echo $(zcat ${RNA_PREFIX}_fpaired.fq.gz |wc -l)/4|bc
    tput setaf 2; echo "------END of  trimming for ${RNA_PREFIX}------"; tput sgr0
    echo
done


# soft mask reads with TETools repeatmodler2 and repeatmasker. 
#    > can struggle if memory is not fast enough and data may need to be stored on compute HPC nodes
#    > This may need to be done more rigourously for larger genomes
~/dfam-tetools.sh
BuildDatabase -name ${DB} ${GENOME}
RepeatModeler -database ${DB} -threads 72 -LTRStruct
RepeatMasker -pa 72 -lib ${DB}-families.fa -xsmall ${GENOME}

# run HiSAT2 or STAR to align (braker3 documentation suggests hisat)
# simplify genome contig names so braker can handle them
IDX=""
GENOME=""
sed -i 's/\s.*$//' ${GENOME}
time hisat2-build ${GENOME} ${IDX}

# EVERYTHING IN ONE ALIGNMENT
NAME=""
FWD_FILES=$(ls -m SRR*fpaired.fq.gz | sed -s 's/ //g')
REV_FILES=$(ls -m SRR*rpaired.fq.gz | sed -s 's/ //g')
time hisat2 -p 32 -q -x ${IDX} -1 ${FWD_FILES} -2 ${REV_FILES} > ${NAME}-hisat2-rnaseq.sam  2> ${NAME}-hisat2-align.err
samtools view -bS -@ 12 ${NAME}-hisat2-rnaseq.sam -o ${NAME}-hisat2-rnaseq.bam
samtools sort -@ 12 ${NAME}-hisat2-rnaseq.bam -o ${NAME}-hisat2-rnaseq_sorted.bam
rm ${NAME}-hisat2-rnaseq.sam; rm ${NAME}-hisat2-rnaseq.bam

# SEPARATE ALIGNMENTS PER FILE
# paired with a list of SRA names (SRRXXXXXX) to align 
for i in $(cat names)
do 
RNASEQ_FWD=${i}_fpaired.fq.gz
RNASEQ_REV=${i}_rpaired.fq.gz
time hisat2 -p 32 -q -x ${IDX} -1 ${RNASEQ_FWD} -2 ${RNASEQ_REV} > ${i}-hisat2-rnaseq.sam  2> ${i}-hisat2-align.err
samtools view -bS ${i}-hisat2-rnaseq.sam -o ${i}-hisat2-rnaseq.bam
samtools sort ${i}-hisat2-rnaseq.bam -o ${i}-hisat2-rnaseq_sorted.bam
rm ${i}-hisat2-rnaseq.sam; rm ${i}-hisat2-rnaseq.bam
done

# unpaired with a list of SRA names (SRRXXXXXX) to align 
for i in $(cat names)
do 
RNASEQ=${i}.fastq.gz
time hisat2 -p 32 -q -x ${IDX} -U ${RNASEQ} > ${i}-hisat2-rnaseq.sam  2> ${i}-hisat2-align.err
samtools view -bS ${i}-hisat2-rnaseq.sam -o ${i}-hisat2-rnaseq.bam
samtools sort ${i}-hisat2-rnaseq.bam -o ${i}-hisat2-rnaseq_sorted.bam
rm ${i}-hisat2-rnaseq.sam; rm ${i}-hisat2-rnaseq.bam
done

# merge bam files
ls unpaired/*sorted.bam >> bam-files
ls paired/*sorted.bam >> bam-files
NAME=$(basename -s .fna ${GENOME}
samtools merge -b bam-files -o ${NAME}-merged.bam --threads 12

# put the aligned bam into BRAKER along with the appropriate orthodb database https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/
T=32
SORTED_BAM="rnaseq_sorted.bam"
HOME=$(ls ~/)
time singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/braker3.sif braker.pl --genome=${GENOME} \
--prot_seq=${HOME}/BRAKER-DB/Alveolata.fa --bam=${SORTED_BAM} --threads=${T} --gff3

# convert the gtf output inro a gff
conda activate agatenv
agat_convert_sp_gxf2gxf.pl -g braker/braker.gtf -o braker/braker.gff
conda deactivate
# compare gffs
gffcompare braker/braker.gff ToxoDB-to-NCBI-ME49-cds-content-check/GCF_000006565.2_TGA4_genomic_d1_f1.0_lifted_annotations-AGAT-phase-fixed_changes_marked.gff -r GCF_000006565.2_TGA4_genomic.gff -o gffcomp-excl-toxo-from-prot-db
gffcompare braker/braker.gff -r ToxoDB-to-NCBI-ME49-cds-content-check/GCF_000006565.2_TGA4_genomic_d1_f1.0_lifted_annotations-AGAT-phase-fixed_changes_marked.gff -o braker-vs-toxodb-excl-toxo-from-protein-db
