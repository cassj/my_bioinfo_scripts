library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
load(filename)

library(ChIPpeakAnno)

#We don't appear to have the current human dataset as a ChIPpeakAnno
#dataset, so we'll have to use the getAnnotation function
library(biomaRt)
ensmart <- useMart("ensembl_mart_43", dataset="hsapiens_gene_ensembl")
TSS.human.NCBI35 <- getAnnotation(ensmart, featureType="TSS")

#need to alter the names so they match ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))
data.annot <- annotatePeakInBatch(rd, AnnotationData=TSS.human.NCBI35)

#and that only gives you the ensembl gene ID, so get extra info:
filters <- c("ensembl_gene_id")
values<-unique(as.character(values(data.annot)[,"feature"]))
attributes <- c("ensembl_gene_id","hgnc_symbol", "description")

more.annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol
more.annot <- more.annot[more.annot[,"hgnc_symbol"]!="",]

#and some have more than one symbol, but the same desc, so ditch them
more.annot <- more.annot[-1*(which(duplicated(more.annot[,c(1,3)]))),]

rownames(more.annot) <- more.annot[,"ensembl_gene_id"]

#and reset the values data.annot with the extra annotation
#and the original scores.
annot <- values(data.annot)
rd <- values(rd)
new.annot <- list()
for(i in 1:length(annot)){
 this.df <- annot[[i]]
 this.rd.df <- rd[[i]]
 these.ids <- this.df[,"feature"]
 this.df<-cbind(this.df,this.rd.df,DataFrame( more.annot[these.ids,2:3]))
 new.annot[[i]] <- DataFrame(this.df)
}

new.annot <- SplitDataFrameList(new.annot, compress=T)
names(new.annot) <- names(annot)

values(data.annot) <- new.annot


newfile<-sub("RangedData", "AnnoRangedData", filename)
save(data.annot, file=newfile)
