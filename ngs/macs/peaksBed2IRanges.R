source("scripts/qw.R")
library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))

skip <- 17
if( length(grep('negative', filename))>0 ){
 skip <- 0
}

data <- read.csv(filename, header=T, sep="\t", skip=skip)

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

outfile=sub('.xls','.RangedData.R',filename)

#And save the result
save(rd, file=outfile)

