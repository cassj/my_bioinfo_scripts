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
#http://www.compbio.group.cam.ac.uk/Resources/Annotation/IlluminaMouseV1.1.txt
#Mouse WG version 1.1 (Mouse-6_v1_1_11234304_A) against Mouse Feb. 2006 (mm8) assembly 
old.annot<-read.csv("lib/IlluminaMouseV1.txt", sep="\t", header=T)

#Our IlluminaIDs correspond to col "Target0" in the annot. Which, means they aren't actually
#Probe IDs, they're transcript IDs. So we should go with the genome position of the transcript.
#The results have been (I assume) averaged over all probes that hit that gene - which is stupid,
#but we don't have access to the raw data so there's nothing we can do about it.
#

#surely illumina *must* have annotation for this, even if it's shite.
# wget http://www.switchtoi.com/pdf/Annotation%20Files/Mouse/Mouse-6_V1.csv.zip

illumina.annot<-read.csv("lib/Mouse-6_V1.csv")

#hmm, doesn't get me any genome positions though, just Accession numbers. Try them in biomaRt


library(biomaRt)
ensmart <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
filters <- 'illumina_mousewg_6_v1'
filters <- 'embl'
values <- as.character(illumina.annot[1:10,"Accession"])
attributes <- qw(embl, ensembl_gene_id, ensembl_transcript_id, start_position, end_position, strand, mgi_symbol)

res <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)




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

             



