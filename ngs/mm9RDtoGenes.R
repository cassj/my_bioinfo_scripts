library(IRanges)
source("scripts/qw.R")

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

data.annot <- as.data.frame(data.annot)
rownames(data.annot) <- as.character(data.annot$names)

#and that only gives you the ensembl gene ID, so get extra info:
library(biomaRt)

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
                
filters <- c("ensembl_gene_id")
values<-unique(as.character(data.annot[,"feature"]))
attributes <- c("ensembl_gene_id","mgi_symbol", "description")

more.annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol
more.annot <- more.annot[more.annot[,"mgi_symbol"]!="",]

#and some have more than one symbol, but the same desc, so ditch them
more.annot <- more.annot[-1*(which(duplicated(more.annot[,c(1,3)]))),]

rownames(more.annot) <- more.annot[,"ensembl_gene_id"]

#add this extra data to the data.annot
data.annot <- data.frame(data.annot, more.annot[data.annot[,"feature"],])

#and add all of this to your original rd
rd<-as.data.frame(rd)
nms<-as.character(rd$names)
rownames(rd)<-nms
data.annot<-data.annot[nms,]   #sort annotation to same order as rd
rd<-cbind(rd, data.annot) 

colnames(rd)<-gsub( "values.","", colnames(rd))
colnames(rd)[17] <- "gene.strand"
colnames(rd)[19] <- "gene.start"
colnames(rd)[20] <- "gene.end"
colnames(rd)[21] <- "is.inside.gene"
colnames(rd)[22] <- "distance.to.gene"
colnames(rd)[23] <- "ensembl.gene.id"
colnames(rd)[24] <- "gene.mgi.symbol"
colnames(rd)[25] <- "gene.description"

values<-qw(Length, Summit, nTags, neg10log10pVal, FoldEnrichment, FDR, gene.strand, gene.start, gene.end, is.inside.gene, distance.to.gene, ensembl.gene.id, gene.mgi.symbol, gene.description) 

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
