#need to make this into a package.
library(IRanges)
library(Biostrings)
library(seqLogo)
library(BSgenome.Mmusculus.UCSC.mm9)


#TODO: How to handle masking?


#get sample regions from a genome
sample.genome <- function(genome,nsamples, nbases, chrs=NULL){
  
  if(is.null(chrs)){chrs <- seqnames(genome)}

  #get chr lengths
  chr.lengths <- seqlengths(genome)[chrs]

  # weight chrs for sampling by seq length
  genome.len <- sum(as.numeric(chr.lengths))
  chr.weights <- chr.lengths/genome.len

  # choose chrs to sample from 
  sample.chrs <- sample(chrs, nsamples, replace=T, prob=chr.weights)
  sample.counts <- table(sample.chrs)

  
  # get the sample sequences
  for(chr in chrs{

    starts <- sample(1:chr.lengths[chr], sample.counts[chr], replace=F)
    ends <- starts+(nbases-1)

    #this drops all rpt-mask masks, but includes regions with undef seqs.
    v <- Views(genome[[chr]],starts, ends)
    #so what do I do with sequences which don't have a 

   }

}

#this assumes Transfac/JASPAR format
read.pfm <- function(file){
  pfm <- as.matrix(read.table(file))
  if (
  rownames(pfm) <- c('A', 'C', 'G', 'T')
  colnames(pfm) <- paste("pos",1:ncol(pfm), sep="")
  return(pfm)
}


#many motifs suffer from small sample biases, add a pseudocount to cope with this
#counts are added by position (count/4 added to each letter)
#can think of it as a bayesian prior if you feel so inclined.
add.pseudocount <- function(pfm, count=1){
  return(pfm + (pseudocount/4))
}

#convert relative frequencies to probabilites
pfm.to.ppm <- function(pfm){
  #just does: pfm[,i]/sum(pfm[,i]) for each col)
  return(prop.table(pfm, margin=2))
}

  
#sequences should be a DNAStringSet
bg.freqs <- function(sequences){

  ab <- c("A","C","G","T")

  tot <- sum(nchar(sequences))
  B <- lapply(ab, function(x){
    sum(vcountPattern(x,sequences))/tot
  })
  names(B) <- ab
  return(B)
}

ppm.to.pwm <- function(ppm, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){
  #careful if we have 0 in the ppm, we'll get Inf
  if(length(which(ppm==0))>0){
    warning("ppms containing zeros will generate Inf values in pwm. Consider pseudocounts")
    }
  pwm <- log2(ppm/bg)
  
  return(pwm)
}


#information content for each col. 
ppm.to.ic <- function(ppm, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)){

  #careful if we have 0 in the ppm, we'll get NaN
   if(length(which(ppm==0))>0){
     warning("ppms containing zeros will generate NaN IC values.  Consider pseudocounts")
    }

  ic <- apply(ppm * ppm.to.pwm(ppm),2,sum)
  return(ic)
}

  


  
