#!/bin/bash

## loading modules
module load umitools
module load seqtk
module load blast
module load samtools
module load R

## Step 1--Download *_CAR_MiSeq fastq files of each sample from GEO

## Step 2 Whitelist generation from read R1 of 5' gene expression data

# generate Sample_CTRL_0h whitelist

umi_tools whitelist --stdin Sample_CTRL_0h_RNA_H25HKBGXC_H5WJHBGXC_H7H3GBGXC_HVJVWBGXB_HVWJ7BGXB_R1.fastq \
                      --plot-prefix=whitelist/CTRL_0h_plot \
                      --set-cell-number=5000 \
                      --error-correct-threshold=1 \
                      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                      --log2stderr > whitelist/CTRL_0h_whitelist.txt

# generate Sample_CTRL_24h whitelist
umi_tools whitelist --stdin Sample_CTRL_24h_RNA_H25HKBGXC_H5WJHBGXC_H7H3GBGXC_HVJVWBGXB_HVWJ7BGXB_R1.fastq \
                      --plot-prefix=whitelist/CTRL_24h_plot \
                      --set-cell-number=15000 \
                      --error-correct-threshold=1 \
                      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                      --log2stderr > whitelist/CTRL_24h_whitelist.txt

# generate Sample_STIM_24h whitelist
umi_tools whitelist --stdin Sample_STIM_24h_RNA_H25HKBGXC_H5WJHBGXC_H7H3GBGXC_HVJVWBGXB_HVWJ7BGXB_R1.fastq \
                      --plot-prefix=whitelist/STIM_24h_plot \
                      --set-cell-number=15000 \
                      --error-correct-threshold=1 \
                      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                      --log2stderr > whitelist/STIM_24h_whitelist.txt

## Step 3-- Extract cell barcode and UMI
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin car_fastq/Sample_CTRL_0h_CAR_MiSeq_D8D5T_D9DB8_JB697_R1.fastq.gz \
                  --stdout extracted_fq/CTRL_0h_CAR_MiSeq_R1_extracted.fastq.gz \
                  --read2-in car_fastq/Sample_CTRL_0h_CAR_MiSeq_D8D5T_D9DB8_JB697_R2.fastq.gz \
                  --read2-out extracted_fq/CTRL_0h_CAR_MiSeq_R2_extracted.fastq.gz \
                  --error-correct-cell \
                  --whitelist=whitelist/CTRL_0h_whitelist.txt


umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin car_fastq/Sample_CTRL_24h_CAR_MiSeq_CKMN5_D9DB8_JB697_R1.fastq.gz \
                  --stdout extracted_fq/CTRL_24h_CAR_MiSeq_R1_extracted.fastq.gz \
                  --read2-in car_fastq/Sample_CTRL_24h_CAR_MiSeq_CKMN5_D9DB8_JB697_R2.fastq.gz \
                  --read2-out extracted_fq/CTRL_24h_CAR_MiSeq_R2_extracted.fastq.gz \
                  --error-correct-cell \
                  --whitelist=whitelist/CTRL_24h_whitelist.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin car_fastq/Sample_STIM_24h_CAR_MiSeq_CKMN5_D9DB8_JB697_R1.fastq.gz \
                  --stdout extracted_fq/STIM_24h_CAR_MiSeq_R1_extracted.fastq.gz \
                  --read2-in car_fastq/Sample_STIM_24h_CAR_MiSeq_CKMN5_D9DB8_JB697_R2.fastq.gz \
                  --read2-out extracted_fq/STIM_24h_CAR_MiSeq_R2_extracted.fastq.gz \
                  --error-correct-cell \
                  --whitelist=whitelist/STIM_24h_whitelist.txt


## Step 4-- Using blastn to align

# convert fastq to fasta
gunzip extracted_fq/*fastq.gz

# creat reference of CAR sequences
makeblastdb -in CARgenome.fa -dbtype nucl -parse_seqids -out car_db

seqtk seq -a extracted_fq/CTRL_0h_CAR_MiSeq_R2_extracted.fastq > extracted_fq/CTRL_0h_CAR_MiSeq_R2_extracted.fasta

# or parallel run multiple samples
nano sample_names.txt
CTRL_0h_CAR_MiSeq
CTRL_24h_CAR_MiSeq
STIM_24h_CAR_MiSeq

cat sample_names.txt | parallel "seqtk seq -a extracted_fq/{}_R1_extracted.fastq > extracted_fq/{}_R1_extracted.fasta"
cat sample_names.txt | parallel "seqtk seq -a extracted_fq/{}_R2_extracted.fastq > extracted_fq/{}_R2_extracted.fasta"

# run alignment using blastn
blastn -query extracted_fq/CTRL_0h_CAR_MiSeq_R2_extracted.fasta -db ref/car_db/car -outfmt 6 > CTRL_0h_CAR_MiSeq_R2_blastn.tsv

# or parallel run multiple samples

cat sample_names.txt | parallel "blastn -query extracted_fq/{}_R2_extracted.fasta -db ref/car_db/car -outfmt 6 > {}_R2_blastn.tsv"

# or submit A swarm task

cat blastn.swarm
blastn -db ref/car_db/car -query extracted_fq/CTRL_0h_CAR_MiSeq_R1_extracted.fasta -outfmt 6 -num_threads 4 -out blastn_outs/CTRL_0h_CAR_MiSeq_R1_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/CTRL_0h_CAR_MiSeq_R2_extracted.fasta -outfmt 6 -num_threads 4 -out blastn_outs/CTRL_0h_CAR_MiSeq_R2_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/CTRL_24h_CAR_MiSeq_R1_extracted.fasta -outfmt 6 -num_threads 4 -out blastn_outs/CTRL_24h_CAR_MiSeq_R1_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/CTRL_24h_CAR_MiSeq_R2_extracted.fasta -outfmt 6 -num_threads 4 -out blastn_outs/CTRL_24h_CAR_MiSeq_R2_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/STIM_24h_CAR_MiSeq_R1_extracted.fasta -outfmt 6 -num_threads 4 -out blastn_outs/STIM_24h_CAR_MiSeq_R1_blastn.tsv
blastn -db ref/car_db/car -query extracted_fq/STIM_24h_CAR_MiSeq_R2_extracted.fasta -outfmt 6 -num_threads 4 -out blastn_outs/STIM_24h_CAR_MiSeq_R2_blastn.tsv

swarm -f ../blastcmd.swarm -t 12 -g 100 --gres=lscratch:200 --time 03-00:00:00 --module blast

## Step 5 --parseblast to finalMatrix
Rscript parseBlast_MiSeq.R CTRL_0h_CAR_MiSeq_R2_blastn.tsv outs/CTRL_0h_CAR_MiSeq

