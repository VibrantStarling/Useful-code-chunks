## get funannotate through singularity
# singularity build funannotate.sif docker://nextgenusfs/funannotate:latest

# get your files in order
GENOME=genome.fna
OUT_DIR=outdir_name
SPECIES=Species_name
STRAIN=strain_name
GENEMARK_GTF=gtf_name
ORTHODB=DB_name

## CLEAN
# funannotate apparently doesn't like spaces in the fastq headers, so if it fails complaining that the headers don't pair properly, 
# try replacing all spaces with _
# this might take a WHILE because these files can be BIG
for i in $(ls *fastq.gz)
do
unpigz -c ${i} | sed -e 's/ /_/g' | pigz > temp.gz
mv temp.gz ${i}
done

# make sure your genome has simplified names. Augustus does not deal well with contig names longer than 16 characters
# this code will remove any trailing characters after the first space of you contig names inplace, so remove the -i flag or make a back up if you need to
sed -i 's/\s.*$//' ${GENOME}

# alternatively, funannotate has a sorting script which will rename and sort contigs by length
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif  funannotate sort -i Spades.genome.cleaned.fa -b scaffold -o Spades.genome.sorted.fa

## MASK
# preferably, soft mask reads with TETools repeatmodler2 and repeatmasker. 
#    > can struggle if memory is not fast enough and data may need to be stored on compute HPC nodes
#    > This may need to be done more rigourously for larger genomes
~/dfam-tetools.sh
BuildDatabase -name ${DB} ${GENOME}
RepeatModeler -database ${DB} -threads 72 -LTRStruct
RepeatMasker -pa 72 -lib ${DB}-families.fa -xsmall ${GENOME}

# otherwise, funannotate can do this very simply with tantan (by default, see https://funannotate.readthedocs.io/en/latest/commands.html for options)
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate mask -i Spades.genome.cleaned.fa --cpus 12 -o masked-${GENOME}

## TRAIN
# replace 'ERR*_1.fastq.gz' and 'ERR*_2.fastq.gz' with whatever pattern covers all the SRA libararies you want to include
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate train -i masked-${GENOME} -o ${OUT_DIR} -l ERR*_1.fastq.gz -r ERR*_2.fastq.gz \
--species ${SPECIES} --strain ${STRAIN} --cpus 48 --max_intronlen 100000 

## PREDICT
# train provides the next bit of code to run. Take that and add any extra options like a bam file and a genemark file if genemark is in your path.
# --repeats2evm is good as evidence modeler takes into account repeat structures, good for most eukaryotes and repetetive sequence I think. Also note organism other, without it, i think it uses some presets for fungal species which you maybe dont want
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate predict -i ${GENOME} -o ${OUT_DIR} -s ${SPECIES} --strain ${STRAIN} --cpus 48 \
--repeats2evm --organism other

# this will be update gene models, while adding alternate isoforms and utr information
singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/funannotate.sif funannotate update -i ${OUT_DIR} -g ${GENOME} --cpus 48 --max_intronlen 100000


# BENEFITS OF FUNANNOTATE over braker3
# produces a stringtie file as default
# produces UTR predictions
# more flexible with the amount of evidence it can take
# is not dependant on dead proprietary software, but can use it (genemark)
# very verbose and guides you through its steps

# CONS OF FUNANNOTATE
# requires time consuming fastq preprocessing to remove spaces from headers
# has a manual model fixing step (I have just deleted the first model if there are two overlapping models of the same length)

