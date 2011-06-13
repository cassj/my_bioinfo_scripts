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

#The chromosomes names used in teh NCBIM37 dataset are 1..Y, MT and
#various NT_***** non-standard chrs. 

                                        #remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}


# NOTE: TSS.mouse.NCBIM37 is actually *gene* start and end positions,
# not individual transcripts.

#get the most recent annotation data from ensembl
tss <- getAnnotation(ensmart, "TSS")
save(tss, file=paste(dirname(filename),"/tss.RData",sep=""))

# They will also pull out anything overlapping.
# Post-process if you only want nearest.
nearest.tss.start <- annotatePeakInBatch(rd,
                                         AnnotationData=tss,
                                         PeakLocForDistance = "middle",    # from the middle of the peak
                                         FeatureLocForDistance = "start",  # to the start of the feature
                                         output = "both",
                                         multiple=TRUE
                                         )

# the overlapping stuff would be exactly the same,so just get nearest 
nearest.tss.end <- annotatePeakInBatch(rd,
                                       AnnotationData=tss,
                                       PeakLocForDistance = "middle",    # from the middle of the peak
                                       FeatureLocForDistance = "end",    # to the end of the feature
                                       output = "nearestStart"
                                       )





#and nearest miRNA
mirnas <- getAnnotation(ensmart, "miRNA")
save(mirnas, file=paste(dirname(filename),"/mirna.RData",sep=""))

nearest.mirna.start <- annotatePeakInBatch(rd,
                                           AnnotationData=mirnas,
                                           PeakLocForDistance = "middle",    # from the middle of the peak
                                           FeatureLocForDistance = "start",  # to the start of the feature
                                           output = "both",
                                           multiple=TRUE
                                         )

nearest.mirna.end <- annotatePeakInBatch(rd,
                                         AnnotationData=mirnas,
                                         PeakLocForDistance = "middle",    # from the middle of the peak
                                         FeatureLocForDistance = "end",    # to the end of the feature
                                         output = "nearestStart"
                                         )








#and nearest exon
exons <- getAnnotation(ensmart, "Exon")
save(exons, file=paste(dirname(filename),"/exons.RData", sep="" ))

nearest.exon.start <- annotatePeakInBatch(rd,
                                           AnnotationData=exons,
                                           PeakLocForDistance = "middle",    # from the middle of the peak
                                           FeatureLocForDistance = "start",  # to the start of the feature
                                           output = "both",
                                           multiple=TRUE
                                         )

nearest.exon.end <- annotatePeakInBatch(rd,
                                         AnnotationData=exons,
                                         PeakLocForDistance = "middle",    # from the middle of the peak
                                         FeatureLocForDistance = "end",    # to the end of the feature
                                         output = "nearestStart"
                                         )



# Save everything
save(nearest.tss.start, file=paste(dirname(filename),"/nearest_tss_start.RData", sep=""))
save(nearest.tss.end, file=paste(dirname(filename),"/nearest_tss_end.RData",sep="" ))
save(nearest.mirna.start,file=paste(dirname(filename),"/nearest_mirna_start.RData", sep="" ))
save(nearest.mirna.end, file=paste(dirname(filename),"/nearest_mirna_end.RData", sep=""))
save(nearest.exon.start, file=paste(dirname(filename), "/nearest_exon_start.RData",sep=""))
save(nearest.exon.end, file=paste(dirname(filename),"/nearest_exon_end.RData", sep=""))


# ok, convert everything to dataframes, add values.
# decide what to do with multiple mappings to peaks at a later stage
nearest.tss.start.df<-as.data.frame(nearest.tss.start)
nearest.tss.end.df<-as.data.frame(nearest.tss.end)
nearest.mirna.start.df<-as.data.frame(nearest.mirna.start)
nearest.mirna.end.df<-as.data.frame(nearest.mirna.end)
nearest.exon.start.df<-as.data.frame(nearest.exon.start)
nearest.exon.end.df<-as.data.frame(nearest.exon.end)


#get the rd rownames - these are chr positions.
rd.df <- as.data.frame(rd)
row.nms <- as.character(rd.df$names)
rownames(rd.df) <- row.nms

#stick nearest start and end together into one data frame
nearest.tss.end.df[,"fromOverlappingOrNearest"] = 'NearestEnd'
nearest.mirna.end.df[,"fromOverlappingOrNearest"] = 'NearestEnd'
nearest.exon.end.df[,"fromOverlappingOrNearest"] = 'NearestEnd'

nearest.tss <- rbind(nearest.tss.start.df, nearest.tss.end.df)
nearest.mirna <- rbind(nearest.mirna.start.df, nearest.mirna.end.df)
nearest.exon <- rbind(nearest.exon.start.df, nearest.exon.end.df)


# map nearest data back to original peak data by chromosome position
req.cols <- c('strand', 'feature', 'start_position', 'end_position', 'insideFeature', 'distancetoFeature', 'shortestDistance', 'fromOverlappingOrNearest')

colnames(nearest.tss) <- paste('tss',colnames(nearest.tss), sep=".")
colnames(nearest.mirna) <- paste('mirna', colnames(nearest.mirna), sep=".")
colnames(nearest.exon) <- paste('exon', colnames(nearest.exon), sep=".")

res.tss <-cbind(rd.df[nearest.tss[,"tss.peak"],],nearest.tss[,paste('tss',req.cols,sep=".")])
res.mirna <-cbind(rd.df[nearest.mirna[,"mirna.peak"],],nearest.mirna[,paste('mirna',req.cols,sep=".")])
res.exon <-cbind(rd.df[nearest.exon[,"exon.peak"],],nearest.exon[,paste('exon',req.cols,sep=".")])



# Get extra info from the ensembl ids for transcripts / mirnas (feature IDs are gene names)
filters <- c("ensembl_gene_id")
values<-unique(c(res.tss[,"tss.feature"], res.mirna[,"mirna.feature"]))
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

res.tss <- cbind(res.tss, annot[res.tss[,"tss.feature"],   c("mgi_symbol", "description")])
res.mirna   <- cbind(res.mirna, annot[res.mirna[,"mirna.feature"],   c("mgi_symbol", "description")])


# Still need to map the exons to gene IDs. Can't filter on exons, so have to get genes and exons and build lookup
filters=''
values<-''
attributes <- c("ensembl_gene_id","ensembl_exon_id")
annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

rownames(annot) <- annot[,"ensembl_exon_id"]
res.exon <- cbind(res.exon,ensembl_gene_id=annot[res.exon[,"exon.feature"],"ensembl_gene_id"])

# then get the extra annot from the gene id
filters <- c("ensembl_gene_id")
values<-unique(res.exon[,"ensembl_gene_id"])
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

res.exon  <- cbind(res.exon, annot[res.exon[,"ensembl_gene_id"],   c("mgi_symbol", "description")])

#save the results.
write.csv(res.tss, file=paste(dirname(filename),"/res_tss.csv", sep=""), row.names=F, quote=F)
write.csv(res.mirna, file=paste(dirname(filename),"/res_mirna.csv", sep=""), row.names=F, quote=F)
write.csv(res.exon, file=paste(dirname(filename),"/res_exon.csv", sep=""), row.names=F, quote=F)



























