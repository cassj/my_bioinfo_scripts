#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]

data <- read.csv(filename, header=T, sep="\t", comment.char = "#")

colnames(data) <- qw(Chr, Start, End, Length, Summit, nTags, neg10log10pVal, FoldEnrichment)
values <-  qw(Length, Summit, nTags, neg10log10pVal, FoldEnrichment)

#sometimes they have FDR
if(ncol(data)==9){
  colnames(data)[9] <- "FDR"
  values<-c(values,"FDR")
}

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
data[,"Name"]<-paste(paste(data[,"Chr"],data[,"Start"],sep=":"),data[,"End"], sep="-")


#note - this fails if chr position is undefined, which can happen
#if we've mapped over from somewhere else, so:
ids <- which(is.na(data[,"Chr"]) | is.na(data[,"Start"]) | is.na(data[,"End"]) )
if(length(ids)>0){
  data <- data[-ids,]
}

rd <- RangedData(ranges = IRanges(
                   start= data$Start,
                   end = data$End,
                   names = as.character(data$Name),
                   ),
                 space = as.character(data$Chr),
                 values = data[,values]
                 )

outfile=sub('.xls','.RangedData.RData',filename)

#And save the result
save(rd, file=outfile)

