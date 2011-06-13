#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
rd <- get(load(filename))

library(ChIPpeakAnno)
library(biomaRt)
ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


#load mouse transcripts
data(TSS.mouse.NCBIM37)

#The chromosomes names used in teh NCBIM37 dataset are 1..Y, MT and
#various NT_***** non-standard chrs. 

                                        #remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}

#chippeakanno discards anything in values, so hang on to them
vals <- as.data.frame(rd)
vals <- vals[,grep("values.", colnames(vals), value=T)]
colnames(vals) <- gsub('values.','',colnames(vals))


# NOTE: TSS.mouse.NCBIM37 is actually *gene* start and end positions,
# not individual transcripts.

#get the most recent annotation data from ensembl
tss <- getAnnotation(ensmart, "TSS")
save(tss, file="tss.RData")

#get the nearest start
nearest.tss <- annotatePeakInBatch(rd,
                                   AnnotationData=tss,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "start",
                                   output = "nearestStart"
                                   )

#and nearest end.
nearest.tes <- annotatePeakInBatch(rd,
                                   AnnotationData=tss,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "end",
                                   output = "nearestStart")
                                   

#and nearest miRNA
mirnas <- getAnnotation(ensmart, "miRNA")
save(mirnas, file="mirna.RData")
nearest.mirna <- annotatePeakInBatch(rd,
                                     AnnotationData=mirnas,
                                     PeakLocForDistance="middle",
                                     FeatureLocForDistance="middle",
                                     output ="nearestStart") 


#and later we'll want to know if it overlaps a coding region?
exons <- getAnnotation(ensmart, "Exon")
save(exons, file="exons.RData")

# annotatePeakInBatch is bastard slow, the IRanges findOverlaps function
# works on RangedData and is much faster
overlapping.exon <- findOverlaps(rd, exons, type="any", select="arbitrary", drop=TRUE)
inds <- is.na(overlapping.exon)
overlapping.exon[inds] <- 0
overlapping.exon[!inds] <- 1


#for reasons I don't understand, annotatePeakInBatch doesn't return your
#data in the order you gave it. so, for example rd["1"][1,] is not
#necessarily nearest.tss["1"][1,]. 
rd.df <- as.data.frame(rd)
ord <- as.character(rd.df$names)

nearest.tss <- as.data.frame(nearest.tss)
rownames(nearest.tss) <- as.character(nearest.tss$peak)
nearest.tss <- nearest.tss[ord,]
nearest.tss <- cbind(nearest.tss,overlapping.exon, type="TSS", vals)

nearest.tes <- as.data.frame(nearest.tes)
rownames(nearest.tes) <- as.character(nearest.tes$peak)
nearest.tes <- nearest.tes[ord,]
nearest.tes <- cbind(nearest.tes, overlapping.exon, type="TES", vals)

nearest.mirna <- as.data.frame(nearest.mirna)
rownames(nearest.mirna) <- as.character(nearest.mirna$peak)
nearest.mirna <- nearest.mirna[ord,]
nearest.mirna <- cbind(nearest.mirna, overlapping.exon, type="miRNA", vals)





#and that only gives you the ensembl gene ID, so get extra info:                
filters <- c("ensembl_gene_id")
values<-unique(c(nearest.tss[,"feature"], nearest.tes[,"feature"]))
attributes <- c("ensembl_gene_id","mgi_symbol", "description")

annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol or a description
no.anno <- intersect(
                     which(annot[,"mgi_symbol"]==""),
                     which(annot[,"description"]==""))
annot <- annot[-1*no.anno,]


# a few have multiple bits of annotation
annot<-cbind(annot, alt.annot="")
dups<-annot[duplicated(annot[,"ensembl_gene_id"]), "ensembl_gene_id"]

#keep the first on and add all the others as alt.annot 
for (d in dups){
  inds <- which(annot[,"ensembl_gene_id"]==d)
  this.alt.annot <- annot[inds[-1], c("mgi_symbol", "description")]
  annot[inds[1],"alt.annot"] <- paste(paste(this.alt.annot[,1], this.alt.annot[,2]), collapse="; ")
}
annot <- annot[!duplicated( annot[,"ensembl_gene_id"] ), ]
rownames(annot) <- annot[,"ensembl_gene_id"]


#add the annotation to your nearest info
nearest.tss <-   cbind(nearest.tss,     annot[nearest.tss[,"feature"],   c("mgi_symbol", "description")])
nearest.tes <-   cbind(nearest.tes,     annot[nearest.tes[,"feature"],   c("mgi_symbol", "description")])
nearest.mirna <- cbind(nearest.mirna,   annot[nearest.mirna[,"feature"], c("mgi_symbol", "description")])


# make a table of nearest overall
nearest.nearest <- nearest.tss
smaller <- which( abs(nearest.nearest$distancetoFeature) > abs(nearest.tes$distancetoFeature) )
if (length(smaller)>0) nearest.nearest[smaller,] <- nearest.tes[smaller,]
smaller <- which( abs(nearest.nearest$distancetoFeature) > abs(nearest.mirna$distancetoFeature) )
if (length(smaller)>0) nearest.nearest[smaller,] <- nearest.mirna[smaller,]

#reorder

#colnms <- qw(space, start, end, width, names, peak, strand, feature, start_position, end_position, insideFeature,distancetoFeature, shortestDistance,fromOverlappingOrNearest,overlapping.exon, type, mgi_symbol, description , Length, Summit, nTags,neg10log10pVal, FoldEnrichment, FDR )
#nearest.tss <- nearest.tss[,colnms]
#nearest.tes <- nearest.tes[,colnms]
#nearest.mirna <- nearest.mirna[,colnms]
#nearest.nearest <- nearest.nearest[,colnms]


#save the results as csv and rd?

write.csv(nearest.tss, file="nearest_tss.csv", row.names=F)
write.csv(nearest.tes, file="nearest_tes.csv", row.names=F)
write.csv(nearest.mirna, file="nearest_mirna.csv", row.names=F)
write.csv(nearest.nearest, file="nearest_nearest.csv", row.names=F)









