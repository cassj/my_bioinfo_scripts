#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ChIPpeakAnno)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]

data.annot <- get(load(filename))
#turn it into a data.frame
data<-as.data.frame(data.annot)
#sort by pval
data <- data[order(data$values.neg10log10pVal, decreasing=T),]

number <- args[2]
if (number =="all") number <-nrow(data)

n<-min(nrow(data), number)

#just take top n
data<-data[1:n,]
#and remake into RD

#Fix any incorrectly named chromosomes:
data$space <- as.character(data$space)
data$space <- gsub('MT', 'M', data$space)

values <- grep('values', colnames(data), value=T)
values <- gsub('values.', '', values)
colnames(data) <-  gsub('values.', '', colnames(data))

data.top <- RangedData(ranges = IRanges(
                   start= data$start,
                   end = data$end,
                   names = as.character(data$names),
                   ),
                 space = as.character(data$space),
                 values = data[,values]
                 )

#fetch the peak sequence.
peaksWithSequences = as.data.frame(getAllPeakSequence(data.top, upstream=0, downstream=0, genome=Mmusculus))

rownames(peaksWithSequences) <- as.character(peaksWithSequences$names)
rownames(data) <- as.character(data$names)

#stick the sequences onto your data.frame
peaksWithSequences <- peaksWithSequences[rownames(data),]
data <- cbind(data, sequence = as.character(peaksWithSequences[,'sequence']))

#make the final rd and save it.
values <- c(values,'sequence')
data <- RangedData(ranges = IRanges(
                     start= data$start,
                     end = data$end,
                     names = as.character(data$names),
                     ),
                   space = as.character(data$space),
                   values = data[,values]
                   )

data.df<-as.data.frame(data)

newfile<-sub("peaks", "peaksWithSeqs", filename)
save(data, file=newfile)

newfile<-sub(".R$",".csv",newfile)
write.csv(data.df, file=newfile)
