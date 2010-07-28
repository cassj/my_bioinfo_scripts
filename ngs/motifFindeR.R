#genetic algorithm approach to finding appropriate PSSMs for a list of ranked sequences

#The problem with existing methods is that they either ignore the whole list
#and the rankings, like MEME, or use the whole list but ignore the rankings, like GADEM,
#or use the rankings but only find consensus sequences, like BCRANK.


library(IRanges)
library(seqLogo)
library(Biostrings)
library(bayesm)


###
# parse data
###

data<-get(load("Macs/run4/NA_peaksWithSeqs.AnnoRangedData.R"))
data<-as.data.frame(data)

#The scoring function expects your sequences to be in rank order
data <- data[order(data$values.neg10log10pVal, decreasing=T),]
#get the sequences you want
seqs <- as.character(data[,"values.sequence"])
names(seqs) <- rownames(data)

#parse data for speedy searching later on.
dnaseqs<-DNAStringSet(seqs)




###
# GenAl functions
###



make.alphabet <- function(seqs){
  alphabet<-table(as.factor(unlist(strsplit( paste(seqs, sep="", collapse="")  ,''))))
  alphabet <-alphabet/sum(alphabet)
  nms <- names(alphabet)
  alphabet <- as.numeric(alphabet)
  names(alphabet) <- nms
  return(alphabet)
}


random.pwm <- function(min.W=5, max.W=30, alphabet=c(A=0.25, C=0.25, G=0.25, T=0.25)){
  
  w <- round(runif(1, min=min.W, max=max.W))
  x <- matrix(nrow=length(alphabet), ncol=w)
  rownames(x) <- names(alphabet)
  
  for (j in 1:ncol(x)){
    #pwms are just probs, not log odds cos we're using them to generate sequences
    x[,j] <- rdirichlet( rep(1,nrow(x)) )
  }
  return(x)
}



pwm.makeSeq<-function(pwm){
  ab <- rownames(pwm)
  gen.seq <- apply(pwm, 2, function(x){ab[which(rmultinom(1,1, x)==1 )]})
  return(paste(gen.seq, collapse=""))
}


#init a population of motifs which at least generate sequences
#found in the data.
init.population <- function(size=100, min.W=5, max.W=30, alphabet=c(A=0.25, C=0.25, G=0.25, T=0.25)){
  return(lapply(1:size, function(x){random.pwm(min.W=min.W, max.W=max.W, alphabet=alphabet)}))
}


#score your PWMs against your dataset.
#on the basis of how many times and at what rank, the sequences they
#generate appear in your list
#p is a pwm matrix, as from pop.init
#dnaseqs is a DNAStringSet
score.pwm <- function(p, dnaseqs, nseq=100){

  n <- length(dnaseqs)
  rank.modifier <- 1-((1:n)/n)
  
  #generate a sample of sequences from the PWM
  gen.seqs <- PDict(unlist(lapply(1:nseq ,function(x){pwm.makeSeq(p)})))

  #see how enriched they are in the population of sequences
  match.mat<-vcountPDict(gen.seqs, dnaseqs)

  #match.mat is a matrix of counts of how many times a search seq (rows) appears
  #in each peak seq (cols)
     
  total.matches <- sum(as.numeric(match.mat))
  if (total.matches == 0){return(0)}
  
  #multiply by the rank modifiers
  mod.match.mat <- match.mat*rank.modifier
  match.scores <- sum(as.numeric(mod.match.mat))
  
  ave.mod.rank <- match.scores / total.matches
  return(ave.mod.rank)   
    
 }


score.population <- function(pop, dnaseqs, nseq){
  #calculate scores
  scores <- sapply(pop, function(x){score.pwm(x, dnaseqs, nseq)} )
  return(scores)
}


#mutate an individual
#the mutation rate is in the range 0..1
mutate <- function(p, rate=0.3){
  n <- ncol(p)
  test<-runif(n)
  mut.cols <- which(test <= rate)
  for(i in mut.cols){
    p[,i] <-  rdirichlet( p[,i] )
  }
  return(p)
}


cross.over <- function(one, two){
  #randomly define chromosomal breakpoints
  bp.one <- sample.int(ncol(one),1)
  bp.two <- sample.int(ncol(two),1)
  new <- cbind(one[,1:bp.one], two[,bp.two:ncol(two)])
  return(new)
}


  
breed <- function(pop, weights, cross.rate=0.1, mutation.rate=0.1){
  #take 1 on the basis of weights
  n <- length(pop)
  choices <- apply(rmultinom(n, 1, weights), 2, function(x){which(x==1)})
  
  #take the high scoring stuff, mutate cols a bit
  #what the hell have I done here?
  new.pop <- lapply(choices, function(x){mutate(pop[[x]], rate=mutation.rate)})
  
  #randomly select individuals to cross.
  n <- length(new.pop)
  test<-runif(n)
  cross <- which(test <= cross.rate)

  for(c in cross){
    new.pop[[c]] <- cross.over(new.pop[[c]], new.pop[[sample.int(length(new.pop),1)]])
  }
  return(new.pop)
}





###
# setup GenAl
###


alphabet <- make.alphabet(seqs)
pop <- init.population(size=30, alphabet=alphabet)


###
# iterate
###

for(i in 1:20){
  scores <- score.population(pop, dnaseqs, 10)
  cat(i,' : ', scores, "\n")
  weights <- scores/sum(scores)
  pop <- breed(pop, weights, cross.rate=0.1, mutation.rate=0.1)
}

final.scores <-  score.population(pop, dnaseqs, 10)  
ranks <- order(final.scores, decreasing=T)


#ok, this kinda converges, but somehow we end up with all the cols being 1,0,0,0 or similar.
#why?
