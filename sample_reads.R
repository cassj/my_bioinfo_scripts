#this should be a package to itself


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
