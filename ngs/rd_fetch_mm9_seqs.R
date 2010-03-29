library(IRanges)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ChIPpeakAnno)

#call like R --vanilla --args filename=\"thing\" number=10000
args<-commandArgs()

#retrieve the filename from the command line.
#Note that we're just eval'ing whatever the user supplies.
eval(parse(text=args[grep('filename', args)]))
data.annot <- get(load(filename))

eval(parse(text=args[grep('number', args)]))

n<-min(nrow(data), number)

#turn it into a data.frame
data<-as.data.frame(data.annot)
#sort by pval
data <- data[order(data$values.neg10log10pVal, decreasing=T),]
#just take top n
data<-data[1:n,]
#and remake into RD


values <- grep('values', colnames(data), value=T)
values <- gsub('values.', '', values)
colnames(data) <-  gsub('values.', '', colnames(data))

data.top <- RangedData(ranges = IRanges(
                   start= data$start,
                   end = data$end,
                   names = data$names,
                   ),
                 space = data$space,
                 values = data[,values]
                 )

#fetch the peak sequence.
peaksWithSequences = as.data.frame(getAllPeakSequence(data.top, upstream=0, downstream=0, genome=Mmusculus))

rownames(peaksWithSequences) <- peaksWithSequences$names
rownames(data) <- data$names

#stick the sequences onto your data.frame
peaksWithSequences <- peaksWithSequences[rownames(data),]
data <- cbind(data, sequence = as.character(peaksWithSequences[,'sequence']))

#make the final rd and save it.
values <- c(values,'sequence')
data <- RangedData(ranges = IRanges(
                     start= data$start,
                     end = data$end,
                     names = data$names,
                     ),
                   space = data$space,
                   values = data[,values]
                   )

newfile<-sub("peaks", "peaksWithSeqs", filename)

save(data, file=newfile)
