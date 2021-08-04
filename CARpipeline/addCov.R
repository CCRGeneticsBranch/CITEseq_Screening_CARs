library(stringr)
library(dplyr)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

readFile = as.character(args[1])

myFile = read.delim(readFile,sep = "\t",header = F)

myFile$V11 = myFile$V10

myFile$V11 = gsub("A|T|C|G","~",myFile$V11 )

write.table(myFile,paste0("final_",readFile) , row.names =F,sep = "\t",quote=F,col.names=F)
