#this is a big slow mess - am still figuring out what I want to do with the motifs
#wil figure out how to make it quick and probly turn it into an R package soon.

library(IRanges)
library(rGADEM)
library(MotIV)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
source("scripts/qw.R")

# call like R --vanilla --args filename=\"thing\" motiffile=\"/space/motifs/MA\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
eval(parse(text=args[grep('motiffile', args)]))
eval(parse(text=args[grep('outfile', args)]))

rd <- get(load(filename))

# load the motif (rest is MA0138.2 )
motif <- as.matrix(read.table(motiffile))
rownames(motif) <- c('A', 'C', 'G', 'T')

# read the sequences from the fasta fle
sequences <- read.DNAStringSet(filename, "fasta")



#find the best hit, and it's score, in each seq
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
  f.hits <- matchPWM(motif, seq, min.score='70%')
  best.f <- 0
  if(length(f.hits)>0){
    for(j in 1:length(f.hits)){
      this <- PWMscoreStartingAt(motif, seq, starting.at=start(f.hits[j]))
      if (this > best.f){
        best.f <- this
        best.f.ind <- j
      }
    }
  }
  
  
                                        # reverse hits
  r.hits <- matchPWM(reverseComplement(motif), seq, min.score='70%')
  best.r <- 0
  if(length(r.hits)>0){
    for(j in 1:length(r.hits)){
      this <- PWMscoreStartingAt(reverseComplement(motif), seq, starting.at=start(r.hits[j]))
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


save(sites, file="Top10KSites.R")
