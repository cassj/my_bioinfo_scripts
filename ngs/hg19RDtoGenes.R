library(IRanges)
source("scripts/qw.R")

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
rd <- get(load(filename))

library(ChIPpeakAnno)

#load mouse transcripts
data(TSS.mouse.NCBIM37)

#remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}


#fetch nearest gene annot with ChIPpeakAnno
library(biomaRt)
ensmart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

TSS.human.GRCh37 <- getAnnotation(ensmart, featureType="TSS")

data.annot <- annotatePeakInBatch(rd, AnnotationData=TSS.human.GRCh37)
data.annot <- as.data.frame(data.annot)
rownames(data.annot) <- as.character(data.annot$names)

#and that only gives you the ensembl gene ID, so get extra info:
filters <- c("ensembl_gene_id")
values<-unique(data.annot[,"feature"])
attributes <- c("ensembl_gene_id","hgnc_symbol", "description")

more.annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol
more.annot <- more.annot[more.annot[,"hgnc_symbol"]!="",]

#and some have more than one symbol, but the same desc, so ditch them
more.annot <- more.annot[-1*(which(duplicated(more.annot[,c(1,3)]))),]
rownames(more.annot) <- more.annot[,"ensembl_gene_id"]
more.annot <- more.annot[,which(colnames(more.annot)!="ensembl_gene_id")]


#add this extra data to the data.annot
ord <- as.character(data.annot[,"feature"])
data.annot <- data.frame(data.annot, more.annot[ord,])

#and add all of this to your original rd
rd<-as.data.frame(rd)
nms<-as.character(rd$names)
rownames(rd)<-nms
data.annot<-data.annot[nms,]   #sort annotation to same order as rd
rd<-cbind(rd, data.annot) 

colnames(rd)<-gsub( "values.","", colnames(rd))

#pull out the values cols
rm.cols<-qw(start, end, names, space, width)
values<- colnames(rd)[!colnames(rd) %in% rm.cols]


rd.annot <- RangedData(ranges = IRanges(
                                   start= rd$start,
                                   end = rd$end,
                                   names = rd$names
                                ),
                       space = rd$space,
                       values = rd[,values]
                      )




newfile<-sub("RangedData", "AnnoRangedData", filename)
save(rd.annot, file=newfile)


























