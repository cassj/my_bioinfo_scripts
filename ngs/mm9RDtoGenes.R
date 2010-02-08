library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
load(filename)

library(ChIPpeakAnno)

#load mouse transcripts
data(TSS.mouse.NCBIM37)

#need to alter the names so they match ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

data.annot <- annotatePeakInBatch(rd, AnnotationData=TSS.mouse.NCBIM37)

#and that only gives you the ensembl gene ID, so get extra info:
library(biomaRt)

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
                
filters <- c("ensembl_gene_id")
values<-unique(as.character(values(data.annot)[,"feature"]))
attributes <- c("ensembl_gene_id","mgi_symbol", "description")

more.annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol
more.annot <- more.annot[more.annot[,"mgi_symbol"]!="",]

#and some have more than one symbol, but the same desc, so ditch them
more.annot <- more.annot[-1*(which(duplicated(more.annot[,c(1,3)]))),]

rownames(more.annot) <- more.annot[,"ensembl_gene_id"]

#and reset the values data.annot
annot <- values(data.annot)
new.annot <- list()
for(i in 1:length(annot)){
 this.df <- annot[[i]]
 these.ids <- this.df[,"feature"]
 this.df<-cbind(this.df,DataFrame( more.annot[these.ids,2:3]))
 new.annot[[i]] <- DataFrame(this.df)
}

new.annot <- SplitDataFrameList(new.annot, compress=T)
names(new.annot) <- names(annot)
values(data.annot) <- new.annot


ids <- as.character(values(data.annot)[,"feature"])
more.annot <- more.annot[ids,]

more.annot <- cbind(as.data.frame(values(data.annot)),more.annot)
values(data.annot) <- more.annot
class(more.annot)


newfile<-sub("RangedData", "AnnoRangedData", filename)
save(data.annot, file=newfile)
