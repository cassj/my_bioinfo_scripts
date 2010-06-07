library(IRanges)
library(Biostrings)
library(seqLogo)
source("scripts/qw.R")


# call like R --vanilla --args filename=\"thing\" motiffile=\"/space/motifs/MA\"
args<-commandArgs()

eval(parse(text=args[grep('filename', args)]))
eval(parse(text=args[grep('motiffile', args)]))


#testing
#filename<-"Macs/run4/top10K.fa"
#motiffile<-"/space/motifs/MA0138.2.pfm"   #REST

motif.name <- gsub('\\.pfm','',gsub(".*\\/",'', motiffile))

outdir <- unlist(strsplit(filename, '\\/'))
outdir <- paste(outdir[1:length(outdir)-1], collapse="/")

# load the motif (rest is MA0138.2 )
# note that this is a frequency matrix and needs to be converted to pwm
pfm <- as.matrix(read.table(motiffile))
rownames(pfm) <- c('A', 'C', 'G', 'T')

#add pseudocount to each col (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/)
pseudocount <- 0.8 
pseudo.pfm <- pfm + (pseudocount/4)

#calculate relative frequencies (just does: pfm[,i]/sum(pfm[,i]) for each col)
ppm <- prop.table(pfm, margin=2)

# read the sequences from the fasta fle
sequences <- read.DNAStringSet(filename, "fasta")

#calculate the background distribution of bases from the sequences
A <- sum(vcountPattern('A',sequences))
C <- sum(vcountPattern('C',sequences))
G <- sum(vcountPattern('G',sequences))
T <- sum(vcountPattern('T',sequences))

tot <- A+C+G+T
A <- A/tot
C <- C/tot
G <- G/tot
T <- T/tot

B <- c(A,C,G,T)
names(B) <- qw(A,C,G,T)

#and use the motif and bg probabilities to calculate the log likelihood
pwm <- log(ppm/B)


#and draw the logo, just for reference.
logo.file <- paste(outdir,"/seqlogo_",motif.name,'.png', sep = "")
bitmap(logo.file)
seqLogo(ppm)
dev.off()

#######################
#Now scan the sequences with the PWM and find the best hit, and it's score, in each seq
sites <- data.frame(location=numeric(length(sequences)),
                    score=numeric(length(sequences)),
                    direction=numeric(length(sequences)),
                    sequence=character(length(sequences)),
                    stringsAsFactors=F
                    )


for (i in 1:length(sequences)){
  if (i %% 100 ==0 ){cat(i,"\n")}
  seq <- sequences[[i]]
  
                                        # forward hits
  f.hits <- matchPWM(pwm, seq, min.score='70%')
  best.f <- 0
  if(length(f.hits)>0){
    for(j in 1:length(f.hits)){
      this <- PWMscoreStartingAt(pwm, seq, starting.at=start(f.hits[j]))
      if (this > best.f){
        best.f <- this
        best.f.ind <- j
      }
    }
  }
  
  
                                        # reverse hits
  r.hits <- matchPWM(reverseComplement(pwm), seq, min.score='70%')
  best.r <- 0
  if(length(r.hits)>0){
    for(j in 1:length(r.hits)){
      this <- PWMscoreStartingAt(reverseComplement(pwm), seq, starting.at=start(r.hits[j]))
      if (this > best.r){
        best.r <- this
        best.r.ind <- j
      }
    }
  }

  if(best.f==0 && best.r==0){next}
  
  if(best.f > best.r){
    sites[i,"location"]  <- start(f.hits[best.f.ind])
    sites[i, 'score']    <- best.f
    sites[i,'direction'] <- 1
    sites[i,'sequence']  <- toString(f.hits[best.f.ind])
  }else{
    sites[i,"location"]  <- start(r.hits[best.r.ind])
    sites[i, 'score']    <- best.r
    sites[i,'direction'] <- -1
    sites[i,'sequence']  <- toString(r.hits[best.r.ind])
  }

}


save(sites, file=paste(outdir, "/",motif.name,"_Top10KSites.R",sep=""))


#######
# Generate QC plots of the enrichment of hits as we go down the list

#calc enrichment as a sliding window
n <- 20
num.motifs <- numeric()
ave.score <- numeric()

for(i in 1:nrow(sites)){
  if(i%%n == 0){
    these <- sites[(i-(n-1)):i,]
    num.motifs <- c(num.motifs,sum(these[,"location"]!=0))
    ave.score <- c(ave.score,sum(these[,"score"])/n)
    
  }
}

perc.motifs <- num.motifs/n

#and plot the results
fl <- paste(outdir,"/perc_motifs_",motif.name,'.png', sep = "")
bitmap(fl)
plot(n*(1:length(perc.motifs)),num.motifs,  main=paste("Percentage of Motifs in a sliding window of \nlength",n,"over ranked ChIPseq sites"), xlab="Rank", ylab="Local Percentage of Motifs")
dev.off()

fl <- paste(outdir,"/ave_score_",motif.name,'.png', sep = "")
bitmap(fl)
plot(n*(1:length(ave.score)),ave.score,  main=paste("Average motif score in a sliding window of \nlength",n,"over ranked ChIPseq sites"), xlab="Rank", ylab="Local Average Score")
dev.off()







