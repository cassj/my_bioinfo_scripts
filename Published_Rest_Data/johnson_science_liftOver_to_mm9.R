library(IRanges)
source('scripts/qw.R')
source('scripts/liftOver.R')

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data<- read.table(filename, sep="\t", header=T)

#rename the cols to standards headers:
colnames(data) <- qw(Chr, Start, End, Length, Summit, nTags, neg10log10pVal, FoldEnrichment)
values <-  qw(Length, Summit, nTags, neg10log10pVal, FoldEnrichment)

#sometimes they have FDR
if(ncol(data)==9){
  colnames(data)[9] <- "FDR"
  values<-c(values,"FDR")
}


to.map<-data[,1:3]
colnames(to.map) <- c("chr","start","end")
#map hg17 to hg19
mapped.hg19 <- liftOver(to.map , chain.file="lib/hg17ToHg19.over.chain")
#the hg19 to mm9
mapped.mm9 <- liftOver(mapped.hg19 , chain.file="lib/hg19ToMm9.over.chain")

#and stick the results back in the dataframe
data[,1:3] <- mapped.mm9[,qw(chr,start,end)] 

#save the results
outfile=sub('xls','mm9.xls',filename)
write.table(data, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
