wd <- setwd( "/home/cassj/work/MillHill/MASH1/agilent_chipchip")

#Positions here have been LiftOver'd to mm9 

data<-read.csv("Mash1CoC_Long_081104_raw_mm9.csv")


#ditch everything except the one with the min p val
#and fill in the missing bits.
test<-data.frame()
vals<-data.frame()
for(i in 1:nrow(data)){
   cat("line",i,"\n")
   line<-data[i,]

   if(line[1]==""){
     line[1:5]<-vals
   }else{
     vals<-line[1:5]
   }
   
  if(line[5]==line[7]){
   test<-rbind(test,line)
 }
}

data<-test
write.table(data, file="just_best_probe.csv", quote=FALSE, row.names=FALSE, sep="\t")
save(data, file="just_best_probe.RData")

#the chr positions we're interested in are the centres of the probes, so:
chr <- gsub(':.+','',data$mm9Position)
start <-as.numeric(gsub('.+:(.+)-.+','\\1',data$mm9Position))
end <- as.numeric(gsub('.+-','',data$mm9Position))

#can't use mean here cos the numbers are too big.
pos <- (end-start)/2
pos <- sapply(pos,round)
pos <- start+pos



#get closest ensembl transcripts
source("../mouseTSS.R")
transcripts <- get.closest.tss(chr,pos)
transcripts<-transcripts$transcripts

save(transcripts, file="transcripts.R")
write.table(transcripts, file="agilent_chipchip_transcripts_mapping.csv", sep="\t", row.names=FALSE, quote=FALSE)


#ok, that gets us the nearest transcripts, now we need to get some
#more annotation about them

refseq.transcripts <- get.closest.refseq(chr,pos)
refseq.transcripts<-refseq.transcripts$transcripts

save(refseq.transcripts, file="refseq_transcripts.R")
write.table(refseq.transcripts, file="agilent_chipchip_refseq_transcripts_mapping.csv", sep="\t", row.names=FALSE, quote=FALSE)

pred.transcripts <- get.closest.predicted.tss(chr,pos)
pred.transcripts<-pred.transcripts$transcripts

save(pred.transcripts, file="pred_transcripts.R")
write.table(pred.transcripts, file="agilent_chipchip_pred_transcripts_mapping.csv", sep="\t", row.names=FALSE, quote=FALSE)




setwd(wd)
