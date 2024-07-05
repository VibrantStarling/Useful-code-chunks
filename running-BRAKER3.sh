# BRAKER3 pipeline

'''
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

RNA_PREFIX="SRR6493555"
GENOME="GCF_000006565.2_TGA4_genomic.fna.gz"
DB=ME49

# trim reads with trimmomatic
for FILE in $(ls rnaseq-dir/*)
do
trimmomatic SE -phred33 -threads 32 ${FILE} ${FILE}_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done

# define the trimmed file objects
RNASEQ_FWD=$(${RNA_PREFIX}+"_1.fastq_qc.fq")
RNASEQ_REV=$(${RNA_PREFIX}+"_2.fastq_qc.fq")

# soft mask reads with TETools repeatmodler2 and repeatmasker. 
#    > can struggle if memory is not fast enough and data may need to be stored on compute HPC nodes
#    > This may need to be done more rigourously for larger genomes
~/dfam-tetools.sh
BuildDatabase -name ${DB} ${GENOME}
RepeatModeler -database ${DB} -threads 72 -LTRStruct
RepeatMasker -pa 72 -lib ${DB}-families.fa -xsmall ${GENOME}

# run HiSAT2 or STAR to align (braker3 documentation suggests hisat)
# paired
hisat2 -p 32 -q -x genome-idx -1 ${RNASEQ_FWD} -2 ${RNASEQ_REV} > ${RNA_PREFIX}-hisat2-paired-aligned-rnaseq.sam  2> hisat2-paired-align.err
# unpaired
hisat2 -p 32 -q -x genome-idx -U ${RNASEQ_FASTQ} > ${RNA_PREFIX}-hisat2-unpaired-aligned-rnaseq.sam  2> hisat2-unpaired-align.err

#convert sam to bam
samtools view -bS rnaseq.sam | samtools sort - rnaseq_sorted.bam
rm rnaseq.sam

# put the aligned bam into BRAKER along with the appropriate orthodb database https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/
singularity exec -B ${PWD}:${PWD},~/ ~/braker3.sif braker.pl --genome=${GENOME} --prot_seq=~/BRAKER-DB/Alveolata.fa --bam=${SORTED_BAM}


