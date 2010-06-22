#need to make this into a package.
library(IRanges)
library(Biostrings)

#Definitions:
#pfm: position frequency matrix - counts of letters.
#ppm: position probability matrix (relative frequencies)
#pwm: position weight matrix - log likelihood (log(model prob/bg prob))


#testing
filename<-"Macs/run4/top10K.fa"
motiffile<-"/space/motifs/MA0138.2.pfm"   #REST


#Virtual class for all matrix representations of motifs
setClass("motifMatrix",
         representation=representation(
           "VIRTUAL",
           data="matrix",
           alphabet="character",
           bg="numeric"
           )
         )

setMethod("initialize","motifMatrix", function(.Object, data, alphabet=NULL, bg=NULL){
   data <- as.matrix(data)
   if(!is.numeric(data))
     stop("Data must be numeric")
   
   #if alphabet not defined, use rownames of data matrix
   if(is.null(alphabet)){
     if(is.null(rownames(data))) stop("No alphabet provided and cannot use matrix rownames")
     alphabet <-as.character(rownames(data))
   }else{        
     if(any(duplicated(alphabet)))
       stop("Duplicate characters in alphabet")
     if(is.null(rownames(data)))
       rownames(data) <- alphabet
     #otherwise, check given alphabet and rownames agree
     if(! all(rownames(data) %in% alphabet))
       stop("Some rows in data do not correspond to an element in alphabet")
     if(! all(alphabet %in% rownames(data)))
       stop("Some elements in alphabet have no corresponding row in data")
   }
   
   #if bg not defined, all alphabet elements equiprobable:
   if(is.null(bg)) {
     bg <- rep(1/nrow(data), nrow(data))
     names(bg) <- alphabet
   }
   if(!sum(bg)==1)
     stop("background probabilties do not sum to 1")
   if(! all(names(bg) %in% alphabet))
     stop("Background probabilites defined for elements not found in alphabet")
   if(! all(alphabet %in% names(bg)))
     stop("Some elements in alphabet have no corresponding background probability")

   #everything valid. build object.
   .Object@data <-data
   .Object@alphabet <- alphabet
   .Object@bg <- bg
   return(.Object)
 })


#TODO getters for data, alphabet, bg.




#Position Frequency Matrix
setClass("pfm",
         contains=c("motifMatrix"),
         representation=representation(
           pseudocount = "numeric"
           ))
setMethod("initialize","pfm",
          function(.Object,data, pseudocount=0, add.pseudocount=TRUE, ...){
            if(add.pseudocount)
              data <- data + (pseudocount/nrow(data))
            callNextMethod(.Object, data, ...)
          }
          )


#Position Probability (relative freq) Matrix
setClass("ppm",
         contains=c("motifMatrix"))
setMethod("initialize","ppm",
          function(.Object, data, ...) {
            if(! all(apply(data,2,sum)==1))
              stop("All columns in a ppm must sum to 1")
            callNextMethod(.Object, data, ...)
          }
          )

#Position Weight Matrix
setClass("pwm",
         contains=c("motifMatrix"))
setMethod("initialize","pwm",function(.Object,...) callNextMethod(.Object,...))




#A motif. With multiple representations
setClassUnion("pfmOrNULL",c("pfm","NULL"))
setClassUnion("ppmOrNULL",c("ppm","NULL"))
setClassUnion("pwmOrNULL",c("pwm","NULL"))
setClassUnion("numericOrNULL",c("numeric","NULL"))

setClass("motif",
         representation=representation(
           pfm = "pfmOrNULL",
           ppm = "ppmOrNULL",
           pwm = "pwmOrNULL",
           ic = "numericOrNULL"
           )
         )


#setters - not for public consumption

setGeneric(".ppm<-",
           function(.Object, value) standardGeneric(".ppm<-"))
setReplaceMethod(".ppm",
                 signature=signature("motif", "ppm"),
                 function(.Object,value) {
                     .Object@ppm <- value
                     .Object
                 })


setGeneric(".pwm<-",
           function(.Object, value) standardGeneric(".pwm<-"))
setReplaceMethod(".pwm",
                 signature=signature("motif", "pwm"),
                 function(.Object,value) {
                     .Object@pwm <- value
                     .Object
                 })

setGeneric(".ic<-",
           function(.Object, value) standardGeneric(".ic<-"))
setReplaceMethod(".ic",
                 signature=signature("motif", "numeric"),
                 function(.Object,value) {
                     .Object@pwm <- value
                     .Object
                 })


#getters
setGeneric("pfm",
           function(.Object) standardGeneric("pfm"))
setMethod('pfm',
          signature=signature(.Object="motif"),
          function(.Object){.Object@pfm}
)

setGeneric("ppm",
           function(.Object) standardGeneric("ppm"))
setMethod('ppm',
          signature=signature(.Object="motif"),
          function(.Object){
            if(!is.null(.Object@ppm))
              return(.Object@ppm)
            if(is.null(.Object@pfm))
              return(NULL)
            data <- prop.table(.Object@pfm@data, margin=2)
            .ppm(.Object) <- new("ppm", data=data)
            
            return(ppm(.Object))
          })


setGeneric("pwm",
           function(.Object) standardGeneric("pwm"))
setMethod('pwm',
          signature=signature(.Object="motif"),
          function(.Object){
            if(!is.null(.Object@pwm))
              return(.Object@pwm)
            if(is.null(ppm(.Object)))
              return(NULL)
            data <- log2(ppm(.Object)/bg(ppm))
            .pwm(.Object)<-new("pwm", data=data, bg=bg(ppm))
            return(pwm(.Object))
          })






#getters
setMethod('pfm', 'motif',
          function(.Object){.Object@pfm})


#convert relative frequencies to probabilites
pfm.to.ppm <- function(pfm){
  #just does: pfm[,i]/sum(pfm[,i]) for each col)
  return(prop.table(pfm, margin=2))
}

  
  
ppm.to.pwm <- function(ppm, bg=uniform.bg(ppm)){
    
  #careful if we have 0 in the ppm, we'll get Inf
  if(length(which(ppm==0))>0){
    warning("ppms containing zeros will generate Inf values in pwm. Consider pseudocounts")
    }
  pwm <- log2(ppm/bg)
  
  return(pwm)
}

#information content (entropy relative to background model) for each column 
ppm.to.ic <- function(ppm, bg=uniform.bg(ppm)){

  #careful if we have 0 in the ppm, we'll get NaN
   if(length(which(ppm==0))>0){
     warning("ppms containing zeros will generate NaN IC values.  Consider pseudocounts")
    }
   
  ic <- ppm * ppm.to.pwm(ppm)
  return(ic)
}






####
# Utilities
####

#this should generate a motif
#read transfac formatted pfms (jaspar are same format)
read.transfac <- read.jaspar <- function(file, pseudocount=0){
  pfm <- as.matrix(read.table(file))
  rownames(pfm) <- c('A', 'C', 'G', 'T')
  pfm <- new("pfm",data=pfm, pseudocount=pseudocount)
  return(new("motif",pfm=pfm))
}

  

#calculate background frequencies from an XStringSet
bg.freqs <- function(sequences){

  ab <- c("A","C","G","T")

  tot <- sum(nchar(sequences))
  B <- lapply(ab, function(x){
    sum(vcountPattern(x,sequences))/tot
  })
  names(B) <- ab
  return(B)
}







#draw sequence logo for a ppm
seqlogo <- function(ppm, bg=uniform.bg(ppm), ic.scale=T){    

  if(!valid.ppm(ppm)) stop("ppm is invalid. Columns should sum to 1")

  alphabet <- rownames(ppm)
  
chars <- c("A", "C", "G", "T")
    letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos <- ncol(pwm)
    if (ic.scale) {
        ylim <- 2
        ylab <- "Information content"
        facs <- pwm2ic(pwm)
    }
    else {
        ylim <- 1
        ylab <- "Probability"
        facs <- rep(1, npos)
    }
    wt <- 1
    x.pos <- 0
    for (j in 1:npos) {
        column <- pwm[, j]
        hts <- 0.95 * column * facs[j]
        letterOrder <- order(hts)
        y.pos <- 0
        for (i in 1:4) {
            letter <- chars[letterOrder[i]]
            ht <- hts[letterOrder[i]]
            if (ht > 0) 
                letters <- addLetter(letters, letter, x.pos, 
                  y.pos, ht, wt)
            y.pos <- y.pos + ht + 0.01
        }
        x.pos <- x.pos + wt
    }
    grid.newpage()
    bottomMargin = ifelse(xaxis, 2 + xfontsize/3.5, 2)
    leftMargin = ifelse(yaxis, 2 + yfontsize/3.5, 2)
    pushViewport(plotViewport(c(bottomMargin, leftMargin, 2, 
        2)))
    pushViewport(dataViewport(0:ncol(pwm), 0:ylim, name = "vp1"))
    grid.polygon(x = unit(letters$x, "native"), y = unit(letters$y, 
        "native"), id = letters$id, gp = gpar(fill = letters$fill, 
        col = "transparent"))
    if (xaxis) {
        grid.xaxis(at = seq(0.5, ncol(pwm) - 0.5), label = 1:ncol(pwm), 
            gp = gpar(fontsize = xfontsize))
        grid.text("Position", y = unit(-3, "lines"), gp = gpar(fontsize = xfontsize))
    }
    if (yaxis) {
        grid.yaxis(gp = gpar(fontsize = yfontsize))
        grid.text(ylab, x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = yfontsize))
    }
    popViewport()
    popViewport()
    par(ask = FALSE)

}
  


  

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
