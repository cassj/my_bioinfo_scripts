data<-read.csv("ChIPSeq.csv", skip=1, na.strings='N.A.')

#> colnames(data)
# [1] "Cluster.ID"                                      
# [2] "Cluster.Overlap.Count"                           
# [3] "Cluster.Location"                                
# [4] "Cluster.Span"                                    
# [5] "Overlap.Location"                                
# [6] "Overlap.Span"                                    
# [7] "Known.Gene.Left.Target"                          
# [8] "Known.Gene.Left.Target.Gene.Name"                
# [9] "Known.Gene.Left.Distance.To.Target"              
#[10] "Known.Gene.Left.Target.Exon"                     
#[11] "Known.Gene.Left.Target.Location"                 
#[12] "Known.Gene.Left.Target.Span"                     
#[13] "Known.Gene.Left.Target.Strand"                   
#[14] "Binding.Site.Location.to.Known.Gene.Left.Target" 
#[15] "Known.Gene.Right.Target"                         
#[16] "Known.Gene.Right.Target.Gene.Name"               
#[17] "Known.Gene.Right.Distance.To.Target"             
#[18] "Known.Gene.Right.Target.Exon"                    
#[19] "Known.Gene.Right.Target.Location"                
#[20] "Known.Gene.Right.Target.Span"                    
#[21] "Known.Gene.Right.Target.Strand"                  
#[22] "Binding.Site.Location.to.Known.Gene.Right.Target"
 
#> dim(data)
#[1] 8418   22


#####
# write a bed file of the locations.
# am assuming the overlap location approximates the binding site
#####

#this seems to confuse liftover
#write(paste('track name="ChIP-seq TC" description="ChIP-seqTC binding"',"\n"), file="seqs.bed")

cluster<-as.character(data$Cluster.Location)
overlap<-as.character(data$Overlap.Location)
chr <- gsub('(.+):.*','\\1', cluster)
start <- gsub('(.+)-.*','\\1', overlap)
end <- gsub('.+-(.+)\\(.+','\\1', overlap)
peak <- gsub('.+\\((.+)\\)','\\1', overlap)
name <- paste('"',as.character(data$Overlap.Location),'"', sep="")
score <- data$Cluster.Overlap.Count

overlaps<-cbind(chr,start,end,name,score)

overlaps.bed<-apply(overlaps, 1, function(x){
                                   x<-paste(x,collapse=" ");
                                   return(x)
                                 }
           )
write.table(overlaps.bed, file="beds/overlap_span.bed", quote=FALSE, col.names=FALSE, row.names=FALSE)


#can't have a single position in a bed file and having the same
#start and end fails in liftover. 
peaks<-cbind(chr,peak,as.numeric(peak)+1,name,score)

peaks.bed<-apply(peaks, 1, function(x){
                                   x<-paste(x,collapse=" ");
                                   return(x)
                                 }
           )
write.table(peaks.bed, file="beds/overlap_peak.bed", quote=FALSE, col.names=FALSE, row.names=FALSE)


#And now use liftOver to move them to mm9 positions.


###
# Make mart?
###

#make a mart from the mm9 genome positions?

#in the long term, this is probaly the way to go for integration of
#TFBS sites and so on, but we can just map the genome positions for now































####
# Original annotation
####

left.targets <- unique(as.character(data$Known.Gene.Left.Target))
right.targets <- unique(as.character(data$Known.Gene.Right.Target))
left.targets <- left.targets[which(left.targets!="No annotation yet")]
right.targets <- right.targets[which(right.targets!="No annotation yet")]

#> length(left.targets)
#[1] 4409
#> length(right.targets)
#[1] 4389


write.table(left.targets, file="left_targets.csv", quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)
write.table(right.targets, file="right_targets.csv", quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)

intersect.targets <- intersect(left.targets,right.targets)
write.table(intersect.targets, file="intersect.targets.csv", quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE)

#> length(intersect.targets)
#[1] 2642


# target is gene symbol, may be 'No annotation yet' 
# Name is refseq, undef for those with no annot


#same thing but only for those which have a refseq
left.data<-data[grep('NM_*', data$Known.Gene.Left.Target.Gene.Name),]
right.data<-data[grep('NM_*', data$Known.Gene.Right.Target.Gene.Name),]

# > dim(left.data)
# [1] 5051   22
# > dim(right.data)
# [1] 4676   22

left.targets.rs <- unique(as.character(left.data$Known.Gene.Left.Target))
right.targets.rs <- unique(as.character(right.data$Known.Gene.Right.Target))
intersect.targets.rs <- intersect(left.targets.rs,right.targets.rs)


#not sure how else to do this:
left.rs <- character(length(left.targets.rs))
right.rs <- character(length(right.targets.rs))
intersect.rs <- character(length(intersect.targets.rs))

for (i in 1:length(left.targets.rs)) {
  left.rs[i] <- as.character(left.data[ which(left.data$Known.Gene.Left.Target == left.targets.rs[i]), 'Known.Gene.Left.Target.Gene.Name'][1])
  if(is.na(left.rs[i])) {stop(paste("erm: ",i)) }
}
                    
for (i in 1:length(right.targets.rs)) {
  right.rs[i] <- as.character(right.data[ which(right.data$Known.Gene.Right.Target == right.targets.rs[i]), 'Known.Gene.Right.Target.Gene.Name'][1])
}
                    
            
for (i in 1:length(intersect.targets.rs)) {
  intersect.rs[i] <- as.character(right.data[ which(right.data$Known.Gene.Right.Target == intersect.targets.rs[i]), 'Known.Gene.Right.Target.Gene.Name'][1])
}
                    

left.targets.rs<-cbind(symbol=left.targets.rs, refseq=left.rs)
right.targets.rs<-cbind(symbol=right.targets.rs, refseq=right.rs)
intersect.targets.rs <- cbind(symbol=intersect.targets.rs, refseq=intersect.rs)

#so about half have refseq ids, again - about 60% intragenic
# > dim(left.targets.rs)
# [1] 2879    2
# > dim(right.targets.rs)
# [1] 2824    2
# > dim(intersect.targets.rs)
# [1] 1699    2


write.table(left.targets.rs, file="left_targets_rs.csv", quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)
write.table(right.targets.rs, file="right_targets_rs.csv", quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)
write.table(intersect.targets.rs, file="intersect_targets_rs.csv", quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)










#I presume this bit actually represents the binding
wid<-data$Overlap.Span
bitmap(file="overlap_span_dens_dist.png", res=300)
plot(density(wid))
dev.off()

#is this just a product of the sequencing read length?
big.wid <- data$Cluster.Span
bitmap(file="cluster_span_dens_dist.png", res=300)
plot(density(big.wid))
dev.off()

bitmap(file="log_cluster_span_dens_dist.png", res=300)
plot(density(log10(big.wid)))
dev.off()



#what about distance to nearest genes?
left.dist <- data$Known.Gene.Left.Distance.To.Target
right.dist <- data$Known.Gene.Right.Distance.To.Target

#remove any that are undef
left.dist <- left.dist[!is.na(left.dist)]
right.dist <- right.dist[!is.na(right.dist)]

dist.union <- c(left.dist, right.dist)

#plot
bitmap(file="distance_to_tss_dens_dist.png", res=300)
opar <- par(mfrow=c(3,1))
plot(density(left.dist))
plot(density(right.dist))
plot(density(dist.union))
dev.off()

#plot logged
bitmap(file="log_distance_to_tss_dens_dist.png", res=300)
opar <- par(mfrow=c(3,1))
plot(density(log10(left.dist)))
plot(density(log10(right.dist)))
plot(density(log10(dist.union)))
dev.off()


# looks like the data has been filtered:

#> max(data$Cluster.Overlap.Count)
#[1] 333
#> min(data$Cluster.Overlap.Count)
#[1] 9










###
# Grab the nearest gene, whether left or right, for each location
##

colnames <- c("Cluster.ID", "Cluster.Overlap.Count", 'Cluster.Location', 'Cluster.Span', 'Overlap.Location', 'Overlap.Span', 'Target', 'Gene.Name', 'Distance.To.Target', 'Target.Exon', 'Target.Location', 'Target.Span', 'Target.Strand', 'Binding.Site.Location.to.Target', 'side')

nearest <- matrix(nrow=nrow(data), ncol=14)
side <- character(nrow(data))
colnames(nearest) <- colnames[1:14]
tmp<-as.matrix(data)

#for (i in 1:nrow(data)){
#  row <- data[i,]
#  nearest[i,1:6] <- as.matrix(row)[1:6]
#
#  if ( sum( is.na( c(row$Known.Gene.Right.Distance.To.Target, row$Known.Gene.Left.Distance.To.Target)) ) == 2){
#    nearest[i,7:14] <- NA
#    side[i] <- 'none'
#    next
#  }
#  
#  if (is.na(row$Known.Gene.Right.Distance.To.Target)){
#    nearest[i,7:14] <- tmp[i,7:14]
#    side[i] <- 'left'
#    next
#  }
#
#  if (is.na(row$Known.Gene.Left.Distance.To.Target)){
#    nearest[i,7:14] <- tmp[i,15:22]
#    side[i] <- 'right'
#    next
#  } 
#
#  #note, if equal we take the left one. prob just intrageneic
#  #with dist 0, so it shouldn't matter
#  if (row$Known.Gene.Left.Distance.To.Target <= row$Known.Gene.Right.Distance.To.Target){
#    nearest[i,7:14] <- tmp[i,7:14]
#    side[i] <- 'left'
#    next
#  }
#  
#  
#  if (row$Known.Gene.Right.Distance.To.Target < row$Known.Gene.Left.Distance.To.Target){
#    nearest[i,7:14] <- tmp[i,15:22]
#    side[i] <- 'right'
#    next
#  }
#
#}
#save(nearest, file="nearest.R")



#this instead
load("nearest.R")

nearest <- cbind(nearest, side=side)
write.table(nearest, file="nearest_genes.csv", quote=FALSE, row.names=FALSE, sep=",")


#ok, and plot some dists of nearest.

bitmap(file="nearest_dist_to_target.png", res=300)
plot(density(as.numeric(nearest[,"Distance.To.Target"])))
dev.off()

#this is wierd. Why are there so many 0? Are they intragenic, or unknown?
foo<-as.numeric(nearest[,"Distance.To.Target"])
plot(density(as.numeric(nearest[,"Distance.To.Target"])))





bitmap(file="nearest_log_dist_to_target.png", res=300)
plot(density(log10(as.numeric(nearest[,"Distance.To.Target"]))))
dev.off()


#and write a file of the 

nearest.gene.symbols<- nearest[,"Target"]

#> length(nearest.gene.symbols)
#[1] 8418
#> length(unique(nearest.gene.symbols))
#[1] 4543

nearest.accessions <- nearest[,"Gene.Name"]
#> length(nearest.accessions)
#[1] 8418
#> length(unique(nearest.accessions))
#[1] 4788

#not quite sure what this means. accessions are diff transcripts maybe,
#which would have the same symbol

#why do we have about half the names duplicated


write.table(nearest.gene.symbols, file="nearest_gene_symbols.csv", quote=FALSE, row.names=FALSE, sep=",")
write.table(nearest.accessions, file="nearest_accessions.csv", quote=FALSE, row.names=FALSE, sep=",")


# everything has a gene symbol but there are lots of duplicates.
# accessions are duplicated too.



#which have big cluster numbers

big.clusters <- nearest[as.numeric(nearest[,"Cluster.Overlap.Count"])> 40,]

#big.clusters[,c("Target","Cluster.Overlap.Count")]


#this is very wierd.

bitmap(file="cluster_overlap_count.png", res=300)
plot(hist(as.numeric(nearest[,"Cluster.Overlap.Count"])))
dev.off()
bitmap(file="log_cluster_overlap_count.png", res=300)
plot(hist(log10(as.numeric(nearest[,"Cluster.Overlap.Count"]))))
dev.off()



#Distribution of peaks per target?

#grab all targets
targets<-nearest[,"Target"]
unique.targets<-unique(nearest[,"Target"])
unique.targets<-sapply(unique.targets, function(x){sum(targets==x)}

write.table(unique.targets, file="peaks_per_target.csv", quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)

#vast majority have only 1. A reasonable number have 2,3,4 etc.
bitmap(file="peaks_per_target_1to10.png", res=300)
hist(unique.targets[unique.targets<10])
dev.off()

#but a few go all the way up to freaky numbers
                       
bitmap(file="peaks_per_target_10plus.png", res=300)
hist(unique.targets[unique.targets>10])
dev.off()


#A530032D15Rik       S100a10        Tdpoz1        Tmem68      AK139008 
#           16            11            13            16            17 
#         Mup1       Tnfrsf8          Rex2      BC108352      AK054210 
#           41            22            90            14            24 
#     BC048648      AK016672       Cyp3a41        Clec4e         Klra1 
#           15            51            15            22            37 
#6430701C03Rik       Clcn4-2         V2r14     MGC117731        Nalp9b 
#           11            37            20            30            14 
#         Pop4     LOC434171         Snrpn A430108E01Rik     LOC436177 
#           23            20            17            39            19 
#2010016B13Rik 4930433N12Rik        Alkbh5     LOC245297 5730507C01Rik 
#           13            24            12            47           135 
#     BC099486          Ighg      AK007163         Gm906          Cts3 
#           34            14            45            16            13 
#       Zfp369      BC048507     LOC544988 1700001E04Rik      AK007159 
#           54            15            52            18            21 
#2610042L04Rik         Sftpd 1700049E17Rik 1700001F09Rik          Mefv 
#           15            19            58            16            11 
#       Dynlt1        Tcp10b      AK138223          Grm4 1700017G21Rik 
#           21            17            15            14           159 
#         Mycs           Ott          Mid1 
#           12            20            21 


                       

#this is wierd - why is there a big gap?
plot(hist(log(unique.targets, 5)))
