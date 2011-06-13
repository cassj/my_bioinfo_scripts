#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
#filename = "MLA_NS_H3K4me3_CMN054_s_2_export_sorted_nodups-W200-G200-islands-summary-FDR1E-3"

data <- read.csv(filename, header=F, sep="\t", comment.char = "#")

colnames(data) <- qw(Chr, Start, End, nTags_ChIP, nTags_Cnt, PValue, FoldEnrichment, FDR)

##calculate neg10 log pvalues
log <- -10*(log10(data[,"PValue"]))

##calculate length of peak
length <- data[,"End"] - data[,"Start"]

##add to dataframe
data_out <- data.frame(data,log, length)
colnames(data_out) <- qw(Chr, Start, End, nTags_ChIP, nTags_Cnt, PValue, FoldEnrichment, FDR, neg10log10pVal, Length)

values <-  qw(Length, nTags_ChIP, nTags_Cnt, neg10log10pVal, FDR, FoldEnrichment)

#sometimes they have FDR
if(ncol(data)==9){
  colnames(data)[9] <- "FDR"
  values<-c(values,"FDR")
}

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
data_out[,"Name"]<-paste(paste(data_out[,"Chr"],data_out[,"Start"],sep=":"),data_out[,"End"], sep="-")


#note - this fails if chr position is undefined, which can happen
#if we've mapped over from somewhere else, so:
ids <- which(is.na(data[,"Chr"]) | is.na(data[,"Start"]) | is.na(data[,"End"]) )
if(length(ids)>0){
  data <- data[-ids,]
}

rd <- RangedData(ranges = IRanges(
                   start= data_out$Start,
                   end = data_out$End,
                   names = as.character(data_out$Name),
                   ),
                 space = as.character(data_out$Chr),
                 values = data_out[,values]
                 )

#outfile=sub('.','.RangedData.RData',filename)

#And save the result
save(rd, file="./SICER.RangedData.RData")

