## 1) Detection of CAR identity from NextSeq, MiSeq and PacBio-seq
These scripts including nextseq.sh, miseq.sh and pacbio.sh are used to process fastq files of CAR binders from three sequencing methods.
Dependencies are not provided but include: CAR sequences genome, linux environment with modules preinstalled.

## 2) Merge all CAR reads count matrix from three sequencing data
Rscript MergeAllseq.R was used to filter cells with less than three reads and combine three matrices obtained from above sequencing methods to get matrix containing cells with unique CAR annotation.
