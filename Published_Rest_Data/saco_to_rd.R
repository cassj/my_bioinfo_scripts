#call like R --vanilla --args filename=\"thing\"
library(IRanges)
source("scripts/qw.R")
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data<- read.table(filename,sep="\t", header=T, as.is=T)

values <-  qw(GSTseq, NumberCollected, ExactMatches, InexactMatches, MatchedSeq)


#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
data[,"Name"]<-paste(paste(data[,"Chr"],data[,"Start"],sep=":"),data[,"End"], sep="-")


#note - this fails if chr position is undefined, which can happen
#if we've mapped over from somewhere else, so:
ids <- which(is.na(data[,"Chr"]) | is.na(data[,"Start"]) | is.na(data[,"End"]) )
if(length(ids)>0){
  data <- data[-ids,]
}

#aaagh, there are duplicate sequences, but I have no idea what that means....

#                GSTseq NumberCollected ExactMatches InexactMatches
#29105 TTGGTCAGGGTGGTCT               1            0              1
#29107 TTGGTCAGGTTCGTCT               1            0              1
#29108 TTGGTCAGGTTGGTCT               1            1              1
#29109 TTGGTCAGGTTGTTCT               1            0              1
#            MatchedSeq     Start Dir  Chr       End                     Name
#29105 TTGGTCAGGTTGGTCT 144742473   r chr3 144742489 chr3:144742473-144742489
#29107 TTGGTCAGGTTGGTCT 144742473   r chr3 144742489 chr3:144742473-144742489
#29108 TTGGTCAGGTTGGTCT 144742473   r chr3 144742489 chr3:144742473-144742489
#29109 TTGGTCAGGTTGGTCT 144742473   r chr3 144742489 chr3:144742473-144742489

#just chuck them for now:
inds<-which(duplicated(data[,"Name"]))
data <- data[-1*inds,]

rd <- RangedData(ranges = IRanges(
                   start= data$Start,
                   end = data$End,
                   names = as.character(data$Name),
                   ),
                 space = as.character(data$Chr),
                 values = data[,values]
                 )





newfile <- 'results/saco_mm9.RangedData.R'
save(rd, file=newfile)




