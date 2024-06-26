# BRAKER3 pipeline

'''
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

curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh
./dfam-tetools.sh

'''

RNA_FASTQ1="SRR6493555_1.fastq.gz"
RNA_FASTQ2="SRR6493555_2.fastq.gz"
GENOME="GCF_000006565.2_TGA4_genomic.fna.gz"
DB=ME49

# trim reads with trimmomatic
trimmomatic SE -phred33 -threads 32 ${RNA_FASTQ_IN}.gz ${OUTPUT}_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

# soft mask reads with TETools repeatmodler2 and repeatmasker. 
#    > can struggle if memory is not fast enough and data may need to be stored on compute HPC nodes
#    > This may need to be done more rigourously for larger genomes
chmod +x ~/dfam-tetools.sh
~/dfam-tetools.sh
BuildDatabase -name ${DB} ${GENOME}
RepeatModeler -database ${DB} -threads 72 -LTRStruct
RepeatMasker -threads 72 -lib ${DB}-families.fa -xsmall ${GENOME}

# Align with HiSAT2 which is included in the BRAKER container (is part of the BRAKER3 pipeline)
singularity exec -B ${PWD}:${PWD} braker3.sif \
braker.pl --genome=${GENOME} --prot_seq=~/BRAKER_DB/alveolata.fa \
--rnaseq_sets_ids=${RNA_FASTQ1},${RNA_FASTQ2} \
--rnaseq_sets_dirs=./


