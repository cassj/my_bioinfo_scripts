#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = FALSE)

filename <- args[1]


qw <- function(...) {
  as.character(sys.call()[-1])
}


liftOver<-function(data, chain.file, ucsc.format=T, chr.col="chr", start.col="start",end.col="end"){

  #data should be a matrix or dataframe with cols for chr, start and end
  #TODO: Or a RangedData / IRange object, 
 
  this<-data.frame(chr=as.character(data[,chr.col]),start=as.numeric(data[,start.col]),end=as.numeric(data[,end.col]), stringsAsFactors=F)

  #Normal counting specifies ranges in 1-based, fully closed form.
  #UCSC specifies ranges in  0-based, half open
  if (!ucsc.format){
      this$start<-this$start-1
  }

  #use the chrstartend as an id so we can map back to the original data
  ids = paste(as.character(data[,chr.col]),data[,start.col],data[,end.col], sep=".")
  this <- cbind(this, ids)
  #If we have duplicate positions, remove them for the liftOver 
  this <- unique(this)
  #we also need to chuck out anything that doesn't have positional data.
  no.chr <- which(is.na(this[,1]))
  no.start <- which(is.na(this[,2]))
  no.end <- which(is.na(this[3]))

  inds <- unique(c(no.chr,no.start,no.end))

  if( length(inds)>0 ){ this <- this[-1*inds,] }
  
  
  ##all this stuff should be a .C() call but I don't have time to make it work just now.
  in.bed <- tempfile()

  #need to watch we don't get scientific notation printed out
  options(scipen=10)
  write.table(this, file=in.bed, sep="\t", row.names=F, col.names=F, quote=F)

  out.bed <- tempfile()
  out.um.bed <- tempfile()
  lo.cmd <- paste("liftOver", in.bed, chain.file, out.bed, out.um.bed)
  system(lo.cmd)

  try(
    new.bed<-read.table(out.bed, sep="\t") 
      ,silent=T)
  try(
    new.um <- read.table(out.um.bed, sep="\t")
      ,silent=T)
  
  #throw away the files
  unlink(c(in.bed, out.bed, out.um.bed))

  if (!exists('new.bed')){stop("No successful mappings")}

  #use the ids as rownames chuck out the column
  #order in the same order as our original data
  #which should stick NAs in for anything that can't be mapped
  rownames(new.bed) <- new.bed[,4]
  new.bed <- new.bed[,-4]
  new.bed <- new.bed[ids,]

  if(!ucsc.format){
   #put the data back to 1-based
   new.bed[,2] <- new.bed[,2]+1
  }

  #replace the new positions in the original dataset
  data[,chr.col] <- new.bed[,1]
  data[,start.col] <- new.bed[,2]
  data[,end.col] <- new.bed[,3]

  #TODO: return some information about the data that won't map
  
  return(data)
  
}



#load both limma and ly results:
limma<-read.csv("filename")
limma<-limma[,-1]
colnames(ky)[1] <-"IlluminaID"
colnames(limma)[1]<-"IlluminaID"

#Annotation from Illumina
lumi.annot<-read.csv("lib/Mouse-6_V1.csv")

#remove anything with many-to-one probe-to-target mapping, there are only a few
dups<-as.character(lumi.annot\$Target[which(duplicated(lumi.annot\$Target))])
lumi.annot<-lumi.annot[-1*which(as.character(lumi.annot\$Target) %in% dups),]
rownames(lumi.annot)<-as.character(lumi.annot\$Target)

#map the TargetID back to a ProbeID using the original Illumina data
ProbeID<-lumi.annot[as.character(limma\$IlluminaID), "ProbeId"]
limma<-cbind(limma, ProbeID)
length(which(is.na(as.character(limma\$ProbeID))))
#[1] 704
limma<-limma[-1*(which(is.na(as.character(limma\$ProbeID)))),]
rownames(limma)<-as.character(limma\$ProbeID)

ProbeID<-lumi.annot[as.character(ky\$IlluminaID),"ProbeId"]
ky<-cbind(ky, ProbeID)
ky<-ky[-1*(which(is.na(as.character(ky\$ProbeID)))),]
rownames(ky)<-as.character(ky\$ProbeID)


#data is sentrix mouse ref 6 according to the paper
#this has old targetIds, that are no longer in use. So, get the probe seqs
#from old reMoat data
old.annot<-read.csv("lib/IlluminaMouseV1.txt", sep="\t", header=T)
rownames(old.annot)<-as.character(old.annot\$ProbeId0)


limma.annot<-cbind(limma, old.annot[as.character(limma\$ProbeID),])
ky.annot<-cbind(ky, old.annot[as.character(ky\$ProbeID),])

#ditch anything for which we have no annotation
remove<-which(is.na(limma.annot\$Chromosome))
if(length(remove)>0){limma.annot<-limma.annot[-remove,]}

remove<-which(is.na(ky.annot\$Chromosome))
if(length(remove)>0){ky.annot<-ky.annot[-remove,]}

#and liftOver the mm8 probe positions
to.map<-limma.annot[,qw(Chromosome, Start, End)]
colnames(to.map) <- c("chr","start","end")
mapped <- liftOver(to.map , chain.file="lib/mm8ToMm9.over.chain")
limma.annot[,qw(Chromosome, Start, End)] <- mapped[,qw(chr,start,end)]

to.map<-ky.annot[,qw(Chromosome, Start, End)]
colnames(to.map) <- c("chr","start","end")
mapped <- liftOver(to.map , chain.file="lib/mm8ToMm9.over.chain")
ky.annot[,qw(Chromosome, Start, End)] <- mapped[,qw(chr,start,end)]

#ditch anything for which we have no region
remove<-which(is.na(limma.annot\$Chromosome))
if(length(remove)>0){limma.annot<-limma.annot[-remove,]}

remove<-which(is.na(ky.annot\$Chromosome))
if(length(remove)>0){ky.annot<-ky.annot[-remove,]}

#we can't use anything mapped to a random chr, so
ids<-grep('random', ky.annot\$Chromosome)
if(length(ids)>0){ky.annot<-ky.annot[(-1*ids),]}

ids<-grep('random', limma.annot\$Chromosome)
if(length(ids)>0){limma.annot<-limma.annot[(-1*ids),]}





#create a RangedData object for the region start and end:
library(IRanges)

#fix colnames
colnames(limma.annot)[c(1,2,5,6,25,26)]<-qw(IlluminaSentrixTargetID, log2FoldChange, pVal, FDR, SpliceJunction, PerfectMatch)
colnames(ky.annot)[c(1,2,3,4,5,23,24)]<-qw(IlluminaSentrixTargetID, log2FoldChange, FoldChange, pVal, FDR, SpliceJunction, PerfectMatch)


#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
rd.limma <- RangedData(ranges = IRanges(
                           start= limma.annot\$Start,
                           end = limma.annot\$End,
                           names = as.character(limma.annot\$IlluminaSentrixTargetID),
                   ),
                 space = as.character(limma.annot\$Chromosome),
                 values = limma.annot[,
                   qw(log2FoldChange, AveExpr,t,pVal, FDR, B, Search_key0, Target0,
                      ProbeId0, Transcript0, Accession0, Symbol0, Definition0, Sequence,
                      Strand, Cytoband, BlastHitType, OtherHits, SpliceJunction, PerfectMatch,
                      Length_match, X.Similarity, Overall.Similarity, GenBanks, Symbols,
                      Description, Comments
                      )]
                 )

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
rd.ky <- RangedData(ranges = IRanges(
                           start= ky.annot\$Start,
                           end = ky.annot\$End,
                           names = as.character(ky.annot\$IlluminaSentrixTargetID)
                   ),
                 space = as.character(ky.annot\$Chromosome),
                 values = ky.annot[,
                   qw(log2FoldChange, FoldChange, pVal, FDR, Search_key0, Target0,
                      ProbeId0, Transcript0, Accession0, Symbol0, Definition0, Sequence,
                      Strand, Cytoband, BlastHitType, OtherHits, SpliceJunction, PerfectMatch,
                      Length_match, X.Similarity, Overall.Similarity, GenBanks, Symbols,
                      Description, Comments
                      )]
                 )





#And save the result
save(rd.limma, file="[% cell_line %]/expression_data/RangedData_Limma.R")
save(rd.ky, file="[% cell_line %]/expression_data/RangedData_KY.R")

write.csv(limma.annot, file="[% cell_line %]/expression_data/limma_annot.csv")
write.csv(ky.annot, file="[% cell_line %]/expression_data/ky_annot.csv")

             



