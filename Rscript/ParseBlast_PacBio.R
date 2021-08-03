library(stringr)
library(dplyr)
library(reshape2)
library(Biobase)

args <- commandArgs(trailingOnly = TRUE)
blast_file <- as.character(args[1])
out_prefix <- as.character(args[2])


r1 = read.delim(blast_file,sep = "\t",header = F)
colnames(r1) = c("ID", "cart", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

r1$read = gsub("(.*ccs.*?)_.*", "\\1", r1$ID)
r1$bc = gsub(".*ccs.*_(.*)_.*", "\\1", r1$ID)
r1$umi = gsub(".*ccs.*_.*_(.*)", "\\1", r1$ID)
r1$type = gsub(".*ccs(.*?)_.*", "\\1", r1$ID)

min_reads=5

topPerRead = aggregate(bitscore ~ read, data = r1, max)

colnames(topPerRead)  =  c("read","bitscore")

merged  = merge(r1, topPerRead, by = c("read","bitscore"))

merged2 = merged[,c("read","bitscore","bc","umi","cart", "type")]
write.table(merged2, paste(out_prefix, "merged2.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

x = merged2[,c("read","bitscore")]

x_uniq = x %>% count(read) %>% filter(n==1)
x_uniq$n =NULL
merged3 = merged2 %>% inner_join(x_uniq, by=c("read"="read"))
write.table(merged3, paste(out_prefix, "merged3.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = T)


counted = as.data.frame(merged3 %>% count(bc, cart))
write.table(counted, paste(out_prefix, "counted.tsv", sep="_"), quote=F, sep="\t", col.names = F,row.names = F)

a1=dcast(counted, bc  ~ cart)
rownames(a1) = a1$bc
a1$bc = NULL 
a1$na_count <- apply(a1, 1, function(x) sum(is.na(x)))

oneCall = a1[a1$na == 13,] 
multicall = a1[a1$na != 13,] 

multicall[is.na(multicall)] <- 0
oneCall[is.na(oneCall)] <- 0
oneCall$na_count = NULL

oneCall$reads_count = rowSums(oneCall)
write.table(oneCall, paste(out_prefix,"oneCall.tsv",sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

multicall$na_count = NULL

write.table(multicall, paste(out_prefix,"multicall.tsv",sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

#o = colnames(oneCall)[apply(oneCall,1,which.max)]

multicall$top1 <- apply(multicall, 1, function(i) sort(i)[ dim(multicall)[2]])


multicall$rest = rowSums(multicall)-(2*multicall$top1)

multicall2 = multicall[multicall$top1 > (2*multicall$rest) & multicall$top1 > 2,]

write.table(multicall2, paste(out_prefix,"multicall2.tsv",sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

multicall2$rest = NULL
multicall2$top1 =NULL
oneCall$reads_count = NULL

final = rbind(multicall2,oneCall)
print(head(final))

finalCallsMatrix  = as.data.frame(cbind(rownames(final), colnames(final)[apply(final,1,which.max)],apply(final,1,max)))
colnames(finalCallsMatrix) <- c("bc","car","read_count")
write.table(finalCallsMatrix, paste(out_prefix, "final.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = F)

finalCallsMatrix <- t(dcast(finalCallsMatrix, bc ~ car))
finalCallsMatrix[is.na(finalCallsMatrix)] <- 0
write.table(finalCallsMatrix, paste(out_prefix, "finalCallReadCountMatrix.tsv", sep="_"), col.names=F, row.names =T,sep = "\t",quote=F)

finalCalls  = cbind(rownames(final), colnames(final)[apply(final,1,which.max)])
colnames(finalCalls) = c("bc","cart")
write.table(finalCalls,paste(out_prefix, "finalCalls", sep="_"),row.names =F,sep = "\t",quote=F)

print(nrow(multicall2))
print(nrow(oneCall))
