library(stringr)
library(dplyr)
library(reshape2)
library(Biobase)
#library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)

r1 = read.delim(as.character(args[1]),sep = "\t",header = F)
out_prefix <- as.character(args[2])
colnames(r1) = c("read", "cart", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

cut_off <- 5

topPerRead = aggregate(bitscore ~ read, data = r1, max)

colnames(topPerRead)  =  c("read","bitscore")

merged2  = merge(r1,topPerRead,by = c("read","bitscore"))
write.table(merged2, paste(out_prefix, "merged2.tsv", sep="_"), quote=F, sep="\t", col.names = T, row.names = F)

x = paste (merged2$read,merged2$bitscore)
merged3 = merged2[isUnique(x),]
write.table(merged3, paste(out_prefix, "merged3.tsv", sep="_"), quote=F, sep="\t", col.names = T, row.names = F)

print(head(merged2))
print(head(merged3))

#merged2 = merge(merged2,bcumi,by = c("read"))


counted = as.data.frame(merged3 %>% count(bc, cart))
write.table(counted, paste(out_prefix, "counted.tsv", sep="_"), quote=F, sep="\t", col.names = T, row.names = F)

print(head(counted))

a1=dcast(counted, bc  ~ cart)
rownames(a1) = a1$bc
a1$bc = NULL 
a1$na_count <- apply(a1, 1, function(x) sum(is.na(x)))

oneCall = a1[a1$na == 13,] 
multicall = a1[a1$na != 13,] 

multicall[is.na(multicall)] <- 0
oneCall[is.na(oneCall)] <- 0
oneCall$na_count = NULL
multicall$na_count = NULL

multicall$top1 <- apply(multicall, 1, function(i) sort(i)[ dim(multicall)[2]])

multicall$rest = rowSums(multicall)-(2*multicall$top1)

multicall2 = multicall[multicall$top1 > (2*multicall$rest) & multicall$top1 >cut_off,]
multicall2$rest = NULL
multicall2$top1 =NULL
write.table(multicall2, paste(out_prefix, "multicall2.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

final = rbind(multicall2,oneCall)
print(head(final))

finalCallsMatrix  = as.data.frame(cbind(rownames(final), colnames(final)[apply(final,1,which.max)],apply(final,1,max)))
colnames(finalCallsMatrix) <- c("bc","car","read_count")
write.table(finalCallsMatrix, paste(out_prefix, "final.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = F)

finalCallsMatrix <- t(dcast(finalCallsMatrix, bc ~ car))
finalCallsMatrix[is.na(finalCallsMatrix)] <- 0
write.table(finalCallsMatrix, paste(out_prefix, "finalMatrix.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

finalCalls  = cbind(rownames(final), colnames(final)[apply(final,1,which.max)])
colnames(finalCalls) = c("bc","cart")
write.table(paste(out_prefix, "finalCall.tsv", sep="_"),row.names =F,sep = "\t",quote=F)

print(nrow(multicall2))
print(nrow(oneCall))

