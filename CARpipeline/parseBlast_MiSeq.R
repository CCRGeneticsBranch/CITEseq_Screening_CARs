library(stringr)
library(dplyr)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

r1_file <- as.character(args[1])
r2_file <- as.character(args[2])
out_prefix <- as.character(args[3])

cut_off <- 5

r1 = read.delim(r1_file,sep = "\t",header = F)
r2 = read.delim(r2_file,sep = "\t",header = F)


r1$read = str_split_fixed(r1$V1,pattern = "_",2)[,1]
r1$id = str_split_fixed(r1$V1,pattern = "_",2)[,2]
r1$bc = str_split_fixed(r1$id,pattern = "_",2)[,1]
r1$umi = str_split_fixed(r1$id,pattern = "_",2)[,2]

r2$read = str_split_fixed(r2$V1,pattern = "_",2)[,1]
r2$id = str_split_fixed(r2$V1,pattern = "_",2)[,2]
r2$bc = str_split_fixed(r2$id,pattern = "_",2)[,1]
r2$umi = str_split_fixed(r2$id,pattern = "_",2)[,2]

#r2$read = r2$V1


merged  = merge(r1,r2,by = c("read","V2"))

merged$overallbit = merged$V12.x + merged$V12.y


topPerRead = aggregate(overallbit ~ read, data = merged, max)

colnames(topPerRead)  =  c("read","overallbit")

merged2  = merge(merged,topPerRead,by = c("read","overallbit"))


merged2 = merged2[,c("read","overallbit","bc.x","umi.x","V2")]

x = merged2[,c("read","overallbit")]

x_uniq = x %>% count(read) %>% filter(n==1)
x_uniq$n=NULL
merged2 = merged2 %>% inner_join(x_uniq, by=c("read"="read"))

#merged2 = merged2[!duplicated(x),]

counted = as.data.frame(merged2 %>% count(bc.x, V2))

write.table(counted, paste(out_prefix, "counted.tsv", sep="_"), quote=F, sep="\t", col.names = F,row.names = F)

a1=dcast(counted, bc.x  ~ V2)
rownames(a1) = a1$bc.x
a1$bc.x = NULL 
a1$na_count <- apply(a1, 1, function(x) sum(is.na(x)))

oneCall = a1[a1$na == 13,] 
multicall = a1[a1$na != 13,] 

multicall[is.na(multicall)] <- 0
oneCall[is.na(oneCall)] <- 0
oneCall$na_count = NULL
multicall$na_count = NULL


#z = colnames(oneCall)[apply(oneCall,1,which.max)]


multicall$top1 <- apply(multicall, 1, function(i) sort(i)[ dim(multicall)[2]])

#sec.large <- as.vector(apply(multicall, 1, order, decreasing=T)[2,])
#print(sec.large)
#multicall$top2 <- sapply(1:length(sec.large), function(i) multicall[i, sec.large[i]])
multicall$rest = rowSums(multicall)-(2*multicall$top1)

multicall = multicall[multicall$top1 > (2*multicall$rest) & multicall$top1 >cut_off,]
multicall$rest = NULL
multicall$top1 =NULL

final = rbind(multicall,oneCall)

finalCallsMatrix  = as.data.frame(cbind(rownames(final), colnames(final)[apply(final,1,which.max)],apply(final,1,max)))
colnames(finalCallsMatrix) <- c("bc","car","read_count")
write.table(finalCallsMatrix, paste(out_prefix, "final.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = F)
finalCallsMatrix <- t(dcast(finalCallsMatrix, bc ~ car))
finalCallsMatrix[is.na(finalCallsMatrix)] <- 0
write.table(finalCallsMatrix, paste(out_prefix, "finalMatrix.tsv", sep="_"), quote=F, sep="\t", col.names = T,row.names = T)

finalCalls  = cbind(rownames(final), colnames(final)[apply(final,1,which.max)])
colnames(finalCalls) = c("bc","cart")
write.table(paste(out_prefix, "finalCall.tsv", sep="_"),row.names =F,sep = "\t",quote=F)

print(nrow(multicall))
print(nrow(oneCall))
