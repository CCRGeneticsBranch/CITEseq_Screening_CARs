library(stringr)
library(dplyr)
library(reshape2)
library(eulerr)
library(RColorBrewer)
library(VennDiagram)

args <- commandArgs(trailingOnly = TRUE)
miseq_file <- as.character(args[1])
nextseq_file <- as.character(args[2])
pacbio_file <- as.character(args[3])
scRNA_file <- as.character(args[4])
out_prefix <- as.character(args[5])

## miseq_file <- "outs/CTRL_0h_CAR_MiSeq_final.tsv"
## nextseq_file <- "outs/CTRL_0h_CAR_NextSeq_final.tsv"
## pacbio_file <- "outs/CTRL_0h_CAR_PacBio_final.tsv"
## scRNA_file <- "outs/CTRL_0h.rna_submtx.tsv"
## out_prefix <- "CTRL_0h_uniq/CTRL_0h"


miseq = read.delim(miseq_file,sep = "\t",header = T)
nextseq = read.delim(nextseq_file,sep = "\t",header = T)
pacbio = read.delim(pacbio_file,sep = "\t",header = T)
scRNA = read.delim(scRNA_file,sep = "\t",header = T)

miseq$car <- gsub("pMH.*_(LH.*)", "\\1", miseq$car)
miseq$car <- gsub("_BBz", "", miseq$car)
miseq$car <- gsub("_B7H3", "", miseq$car)
miseq$car <- gsub("_HL", "-HL", miseq$car)
miseq$car <- gsub("_LH", "-LH", miseq$car)

nextseq$car <- gsub("pMH.*_(LH.*)", "\\1", nextseq$car)
nextseq$car <- gsub("_BBz", "", nextseq$car)
nextseq$car <- gsub("_B7H3", "", nextseq$car)
nextseq$car <- gsub("_HL", "-HL", nextseq$car)
nextseq$car <- gsub("_LH", "-LH", nextseq$car)

#the intersect should be 14
print(length(intersect(intersect(miseq$car,nextseq$car),pacbio$car)))

####################################################
## filter nextseq data supported with only 1 and 2 read exist
print(nrow(nextseq[which(nextseq$read_count==1),]))
print(nrow(nextseq[which(nextseq$read_count>=3),]))

nextseq <- nextseq[which(nextseq$read_count>=3),]

## filter miseq data supported with only 1 and 2 read exist
print(nrow(miseq[which(miseq$read_count==1),]))
print(nrow(miseq[which(miseq$read_count>=3),]))

miseq <- miseq[which(miseq$read_count>=3),]

print(nrow(pacbio[which(pacbio$read_count!=1),]))

## get union_bc among above three matrix
## unique_bc <- unique(c(miseq$bc,nextseq$bc,pacbio$bc))

union_bc <- data.frame("bc"=sort(union(miseq$bc,union(nextseq$bc,pacbio$bc))))
merged <- union_bc %>% dplyr::left_join(miseq, by=c("bc"="bc"))
merged <- merged %>% dplyr::left_join(nextseq, by=c("bc"="bc"))
merged <- merged %>% dplyr::left_join(pacbio, by=c("bc"="bc"))
colnames(merged) <- c("bc","car.miseq", "read_count.miseq","car.nextseq", "read_count.nextseq","car.pacbio", "read_count.pacbio")

## filter nextseq data supported by only 1 read exist
#nrow(merged[which(merged$read_count.nextseq==1),])
#nextseq_filter <- which(merged$read_count.nextseq==1)
#merged$car.nextseq[nextseq_filter] <- NA
#merged$read_count.nextseq[nextseq_filter] <- NA

merged$consensus <- apply(merged %>% dplyr::select(car.miseq, car.nextseq, car.pacbio),1,function(x) names(sort(table(x),decreasing=T))[1])
merged$consensus_count <- apply(merged %>% dplyr::select(car.miseq, car.nextseq, car.pacbio),1,function(x) sort(table(x),decreasing=T)[1])
merged$consensus2nd <- apply(merged %>% dplyr::select(car.miseq, car.nextseq, car.pacbio),1,function(x) names(sort(table(x),decreasing=T))[2])
merged$consensus2nd_count <- apply(merged %>% dplyr::select(car.miseq, car.nextseq, car.pacbio),1,function(x) sort(table(x),decreasing=T)[2])

merged$consensus_read_count <- 0
merged$consensus_read_count[which(merged$car.miseq==merged$consensus)] <- merged$consensus_read_count[which(merged$car.miseq==merged$consensus)] + merged$read_count.miseq[which(merged$car.miseq==merged$consensus)]
merged$consensus_read_count[which(merged$car.nextseq==merged$consensus)] <- merged$consensus_read_count[which(merged$car.nextseq==merged$consensus)] + merged$read_count.nextseq[which(merged$car.nextseq==merged$consensus)]
merged$consensus_read_count[which(merged$car.pacbio==merged$consensus)] <- merged$consensus_read_count[which(merged$car.pacbio==merged$consensus)] + merged$read_count.pacbio[which(merged$car.pacbio==merged$consensus)]

merged$consensus2nd_read_count <- 0
merged$consensus2nd_read_count[which(merged$car.miseq==merged$consensus2nd)] <- merged$consensus2nd_read_count[which(merged$car.miseq==merged$consensus2nd)] + merged$read_count.miseq[which(merged$car.miseq==merged$consensus2nd)]
merged$consensus2nd_read_count[which(merged$car.nextseq==merged$consensus2nd)] <- merged$consensus2nd_read_count[which(merged$car.nextseq==merged$consensus2nd)] + merged$read_count.nextseq[which(merged$car.nextseq==merged$consensus2nd)]
merged$consensus2nd_read_count[which(merged$car.pacbio==merged$consensus2nd)] <- merged$consensus2nd_read_count[which(merged$car.pacbio==merged$consensus2nd)] + merged$read_count.pacbio[which(merged$car.pacbio==merged$consensus2nd)]

merged$consensus_final <- NA
#if consensus > 1(take majority as real annotation) or only one called
consensus1st <- which(merged$consensus_count>1 | is.na(merged$consensus2nd))

## check if there are mis-assigned but can ignore because of using sort decreasing=F
consensus2nd <- which(merged$consensus2nd_count>1)
merged$consensus_final[consensus1st] <- merged$consensus[consensus1st]
merged$consensus_final[consensus2nd] <- merged$consensus2nd[consensus2nd]

## check how many cells can not be annotated
print(head(merged[is.na(merged$consensus_final),]))
print(nrow(merged[is.na(merged$consensus_final),]))  ##16
c=merged[is.na(merged$consensus_final),]
write.table(c, paste(out_prefix, "NA_plot_v5.tsv", sep="_"), col.names = T, row.names = F, sep="\t", quote=F)

## When there are mismatch between two methods, pick more reads supported identity as car annotation
read_count_cutoff <- 2
consensus1st <- which(merged$consensus_count==merged$consensus2nd_count & merged$consensus_read_count > 2*merged$consensus2nd_read_count &  merged$consensus_read_count > read_count_cutoff)
consensus2nd <- which(merged$consensus_count==merged$consensus2nd_count & merged$consensus2nd_read_count > 2*merged$consensus_read_count & merged$consensus2nd_read_count > read_count_cutoff)

## consensus1st <- which(merged$consensus_count==merged$consensus2nd_count & merged$consensus_read_count > read_count_cutoff &  merged$consensus2nd_read_count < read_count_cutoff)
## consensus2nd <- which(merged$consensus_count==merged$consensus2nd_count & merged$consensus_read_count < read_count_cutoff &  merged$consensus2nd_read_count > read_count_cutoff)

merged$consensus_final[consensus1st] <- merged$consensus[consensus1st]
merged$consensus_final[consensus2nd] <- merged$consensus2nd[consensus2nd]

##check how many cells can not be annotated
print(head(merged[is.na(merged$consensus_final),]))
print(nrow(merged[is.na(merged$consensus_final),]))   ###5

##filter unmatched annotation in 3 methods data
merged$consensus_final[which(!is.na(merged$car.miseq) & !is.na(merged$car.nextseq) & !is.na(merged$car.pacbio) & merged$consensus_count==1 & merged$consensus2nd_count==1)] <- NA

# get the read count
merged$consensus_read_count_final <- 0
merged$consensus_read_count_final[which(merged$consensus_final == merged$consensus)] <- merged$consensus_read_count[which(merged$consensus_final == merged$consensus)]
merged$consensus_read_count_final[which(merged$consensus_final == merged$consensus2nd)] <- merged$consensus2nd_read_count[which(merged$consensus_final == merged$consensus2nd)]
merged$consensus_labels <- apply(merged %>% dplyr::select(car.miseq, car.nextseq, car.pacbio,consensus_final),1,function(x) paste(unlist(gsub("car\\.","", c(names(x)[which(as.character(x)==as.character(x)[4])]))),collapse = "&"))
merged$consensus_labels <- gsub("&consensus_final","",merged$consensus_labels)
merged$consensus <- NULL
merged$consensus_count <- NULL
merged$consensus_read_count <- NULL
merged$consensus2nd <- NULL
merged$consensus2nd_count <- NULL
merged$consensus2nd_read_count <- NULL

## filter only 1 read and only exist in nextseq data
## nextseq_filter <- which(merged$consensus_read_count_final==1 & merged$consensus_labels == "nextseq")
## merged$consensus_final[nextseq_filter] <- NA

merged <- merged[which(!is.na(merged$consensus_final)),]

write.table(merged, paste(out_prefix, "merged.tsv", sep="_"), col.names = T, row.names = F, sep="\t", quote=F)

## filter in scRNA library
scRNAbc <- data.frame("bc"=sort(colnames(scRNA)))
scRNAbc$scRNA <- 'Y'
merged <- merged %>% dplyr::left_join(scRNAbc, by=c("bc"="bc"))
merged$scRNA[which(is.na(merged$scRNA))] <- "N"

print(table(merged$scRNA))

merged$scRNA <- as.factor(merged$scRNA)
write.table(table(merged %>% dplyr::select(consensus_labels,scRNA)), paste(out_prefix, "merged_summary.tsv", sep="_"), col.names = NA, row.names = T, sep="\t", quote=F)

#merged$consensus_labels[which(merged$scRNA=="Y")] <- paste(merged$consensus_labels[which(merged$scRNA=="Y")], "scRNA",sep="&")
merged_summary <- merged %>% count(consensus_labels)
write.table(merged, paste(out_prefix, "merged_bc.tsv", sep="_"), col.names = T, row.names = F, sep="\t", quote=F)

## plot three data by venn
venn_data <- merged_summary$n
names(venn_data) <- merged_summary$consensus_labels
vd <- euler(venn_data)
eulerr_options(edges = list(col = "grey"), fontsize = 24)
pdf(paste(out_prefix, "venn.pdf", sep="_"))
plot(vd, quantities=T, edges=c("grey"), fills=c("lightblue", "pink", "green", "grey"))
dev.off()

#head(merged[which(!is.na(merged$car.miseq) & !is.na(merged$car.nextseq) & !is.na(merged$car.pacbio) & merged$consensus_count==1 & merged$consensus2nd_count==1),])
merged <- merged %>% dplyr::filter(scRNA=="Y")
matrix <- t(dcast(merged, bc ~ consensus_final, value.var="consensus_read_count_final"))
matrix[is.na(matrix)] <- 0
write.table(matrix, paste(out_prefix, "finalMatrix.tsv", sep="_"), quote=F, sep="\t", col.names = F,row.names = T)

## plot filtered in scRNA data
merged_summary <- merged %>% count(consensus_labels)
write.table(merged, paste(out_prefix, "merged_bc_filtered.tsv", sep="_"), col.names = T, row.names = F, sep="\t", quote=F)

venn_data <- merged_summary$n
names(venn_data) <- merged_summary$consensus_labels
vd <- euler(venn_data)
eulerr_options(edges = list(col = "grey"), fontsize = 24)
pdf(paste(out_prefix, "venn_filtered.pdf", sep="_"))
plot(vd, quantities=T, edges=c("grey"), fills=c("lightblue", "pink", "green", "grey"))
dev.off()


## check pacbio only with 1 read
pacbio_only = merged[which(merged$consensus_labels=="pacbio"),] ##296
nrow(pacbio_only)
nrow(pacbio_only[which(pacbio_only$read_count.pacbio == 1),]) ##151

## list counts per each car
print(nrow(merged))
print(table(merged$consensus_final))
