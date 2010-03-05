source("scripts/qw.R")
source("scripts/liftOver.R")


#load both limma and ly results:
ky<-read.csv("expression_data/ky_results.csv")
limma <- read.csv("expression_data/limma_results.csv")
limma<-limma[,-1]
colnames(ky)[1] <-"IlluminaID"
colnames(limma)[1]<-"IlluminaID"


#data is sentrix mouse ref 6 according to the paper
#this has old targetIds, that are no longer in use. So, get the probe seqs
#from old reMoat data
old.annot<-read.csv("lib/IlluminaMouseV1.txt", sep="\t", header=T)

#a few are duplicated
#old.annot[duplicated(as.character(old.annot[,2])),"Symbols"]
# [1] Txn1 Slc25a33            Txn1 Slc25a33            Eef1a1                  
# [4] Eef1a1                   Eef1a1                   LOC380687 BC096042 Gapdh
# [7] LOC380687 BC096042 Gapdh LOC380687 BC096042 Gapdh Ubc                     
#[10] Ubc                      Ubc                      AK013903 Rps9           
#[13] AK013903 Rps9            AK013903 Rps9            Tubb2b                  
#[16] Actb                     Actb              

#but very few, so let's just remove the dups and go with the probe locations:
old.annot<-old.annot[!duplicated(old.annot[,2]),]
rownames(old.annot)<-old.annot[,2]
limma.annot<-cbind(limma, old.annot[as.character(limma$IlluminaID),])
ky.annot<-cbind(ky, old.annot[as.character(ky$IlluminaID),])

#ditch anything for which we have no annotation
remove<-which(is.na(limma.annot$Chromosome))
limma.annot<-limma.annot[-remove,]

remove<-which(is.na(ky.annot$Chromosome))
ky.annot<-ky.annot[-remove,]

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
remove<-which(is.na(limma.annot$Chromosome))
limma.annot<-limma.annot[-remove,]

remove<-which(is.na(ky.annot$Chromosome))
ky.annot<-ky.annot[-remove,]


#create a RangedData object for the region start and end:
library(IRanges)

#fix colnames
colnames(limma.annot)[c(1,2,5,6,24,25)]<-qw(IlluminaSentrixTargetID, log2FoldChange, pVal, FDR, SpliceJunction, PerfectMatch)
colnames(ky.annot)[c(1,2,3,4,5,22,23)]<-qw(IlluminaSentrixTargetID, log2FoldChange, FoldChange, pVal, FDR, SpliceJunction, PerfectMatch)


#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
rd.limma <- RangedData(ranges = IRanges(
        	           start= limma.annot$Start,
                	   end = limma.annot$End,
	                   names = as.character(limma.annot$IlluminaSentrixTargetID),
                   ),
                 space = as.character(limma.annot$Chromosome),
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
                      start= ky.annot$Start,
                      end = ky.annot$End,
                      names = as.character(ky.annot$IlluminaSentrixTargetID),

                   ),
                    space = as.character(ky.annot$Chromosome),
                    values = ky.annot[,
                   qw(log2FoldChange, FoldChange, pVal, FDR, Search_key0, Target0,  
                      ProbeId0, Transcript0, Accession0, Symbol0, Definition0, Sequence,
                      Strand, Cytoband, BlastHitType, OtherHits, SpliceJunction, PerfectMatch,
                      Length_match, X.Similarity, Overall.Similarity, GenBanks, Symbols,
                      Description, Comments
                      )]
                 )



#And save the result
save(rd.limma, file="expression_data/RangedData_Limma.R")
save(rd.ky, file="expression_data/RangedData_KY.R")

             


