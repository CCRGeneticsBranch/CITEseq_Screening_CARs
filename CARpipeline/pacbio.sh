#!/bin/bash

## loading modules
module load picard
module load umitools
module load seqtk
module load blast
module load samtools
module load R

## Step 1-- Download *_CAR_PacBio fastq files of each sample from GEO

## Step 2-- convert fastq files to sam files

java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar FastqToSam  F1=Sample_CTRL_0h_CAR_PacBio.ccs.fastq O=CTRL_0h_CAR_PacBio_unaligned.sam SM=CTRL_0h_CAR_PacBio

## Step3-- run toSearch.py to organize sequence
python toSearch.py CTRL_0h_CAR_PacBio_unaligned.sam > CTRL_0h_CAR_PacBio_trimmed.sam

# or Split to speed task:
cp CTRL_0h_CAR_PacBio_unaligned.sam CTRL_0h_CAR_PacBio/
split -l 1000 CTRL_0h_CAR_PacBio_unaligned.sam
for i in `ls *`;do echo "python CARpipe/toSearch.py $i >> CTRL_0h_CAR_PacBio_trimmed.sam" >> CTRL_0h_CAR_PacBio.swarm;done
swarm -f CTRL_0h_CAR_PacBio.swarm -g 100 --gres=lscratch:200 --time 03-00:00:00 --module python --logdir swarm_log

## Step4-- convert above sam file to fastq
## add quality score to trimmed.sam file
Rscript CARpipe/addCov.R CTRL_0h_CAR_PacBio_trimmed.sam

java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar SamToFastq I=final_CTRL_0h_CAR_PacBio_trimmed.sam FASTQ=CTRL_0h_CAR_PacBio_trimmed.fastq

## Step 5-- extract cell barcode and UMI
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin trimmed_fastq/S921/CTRL_0h_CAR_PacBio_trimmed.fastq \
                  --stdout extracted_fq/CTRL_0h_CAR_PacBio_extracted.fastq \
                  --error-correct-cell \
                  --whitelist=whitelist/CTRL_0h_whitelist.txt

## Step 6--blastn to align

## make database for car references
makeblastdb -in CARgenome.fa -dbtype nucl -parse_seqids -out car_db/car

seqtk seq -a CTRL_0h_CAR_PacBio_extracted.fastq > CTRL_0h_CAR_PacBio_extracted.fasta

or submit swarm task
cat blastn.swarm
blastn -db ref/car_db/car -query extracted_fq/CTRL_0h_CAR_PacBio_extracted.fasta -outfmt 6 -num_threads 6 -out blastn_outs/CTRL_0h_CAR_PacBio_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/CTRL_24h_CAR_PacBio_extracted.fasta -outfmt 6 -num_threads 6 -out blastn_outs/CTRL_24h_CAR_PacBio_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/STIM_24h_CAR_PacBio_extracted.fasta -outfmt 6 -num_threads 6 -out blastn_outs/STIM_24h_CAR_PacBio_blastn.tsv

swarm -f blastn.swarm -t 12 -g 100 --gres=lscratch:200 --time 03-00:00:00 --module blast --logdir swarm_log

## Step 7--parseblastn to finalMatrix
Rscript parseBlast_PacBio.R blastn_outs/CTRL_0h_CAR_PacBio_blastn.tsv outs/CTRL_0h_CAR_PacBio_blastn.tsv
