#need to make this into a package.
#library(IRanges)
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
           background="numeric"
           )
         )

setMethod("initialize","motifMatrix", function(.Object, data, alphabet=NULL, background=NULL){
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
   if(is.null(background)) {
     background <- rep(1/nrow(data), nrow(data))
     names(background) <- alphabet
   }
   if(!sum(background)==1)
     stop("background probabilties do not sum to 1")
   if(! all(names(background) %in% alphabet))
     stop("Background probabilites defined for elements not found in alphabet")
   if(! all(alphabet %in% names(background)))
     stop("Some elements in alphabet have no corresponding background probability")

   #everything valid. build object.
   .Object@data <-data
   .Object@alphabet <- alphabet
   .Object@background <- background
   return(.Object)
 })


#getters
setGeneric("data",
           function(.Object) standardGeneric("data"))
setMethod("data",
          signature=signature("motifMatrix"),
          function(.Object) {
            .Object@data
          })

setGeneric("alphabet",
           function(.Object) standardGeneric("alphabet"))
setMethod("alphabet",
          signature=signature("motifMatrix"),
          function(.Object) {
            .Object@alphabet
          })

setGeneric("background",
           function(.Object) standardGeneric("background"))
setMethod("background",
          signature=signature("motifMatrix"),
          function(.Object) {
            .Object@background
          })

setGeneric("bg",
           function(.Object) standardGeneric("bg"))
setMethod("bg",
          signature=signature("motifMatrix"),
          function(.Object) {
            .Object@background
          })






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

setGeneric("pseudocount",
           function(.Object) standardGeneric("pseudocount"))
setMethod("pseudocount",
          signature=signature("motifMatrix"),
          function(.Object) {
            .Object@data
          })





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
setClassUnion("characterOrNULL",c("character","NULL"))

setClass("motif",
         representation=representation(
           pfm = "pfmOrNULL",
           ppm = "ppmOrNULL",
           pwm = "pwmOrNULL",
           ic = "numericOrNULL",
           alphabet="characterOrNULL"
           )
         )


setMethod("initialize","motif",
          function(.Object, pfm=NULL, ppm=NULL){
            
            .Object@pfm <- pfm
            .Object@ppm <- ppm
            
            if(!is.null(.Object@pfm)){ #got a pfm
              if(!is.null(.Object@ppm)){ #also got a ppm?
                if(! (all(alphabet(pfm) %in% alphabet(ppm)) && all(alphabet(ppm) %in% alphabet(pfm) )))
                  stop("mismatching alphabets in supplied pfm and ppm")
                warning("Using supplied pfm and ppm. Might be better to provide only pfm, from which ppm can be derived.")
                .Object@alphabet <- alphabet(pfm)
              }
            }else{ 
              if(!is.null(.Object@ppm)){ #just got a ppm
                .Object@alphabet <- alphabet(ppm)
              } else{
                stop("Please provide either a pfm or pwm object") 
              }
            }
            .Object
          }
          )

#setters - not for public consumption

setGeneric(".ppm<-",
           function(.Object, value) standardGeneric(".ppm<-"))
setReplaceMethod(".ppm",
                 signature=signature("motif", "ppm"),
                 function(.Object,value) {
                     if(!is.null(alphabet(.Object))){

                     }
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
                   .Object@ic <- value
                   .Object
                 })


#getters


setGeneric("alphabet",
           function(.Object) standardGeneric("alphabet"))
setMethod('alphabet',
          signature=signature(.Object="motif"),
          function(.Object){.Object@alphabet}
)

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
            data <- log2(data(ppm(.Object))/bg(ppm(.Object)))
            .pwm(.Object)<-new("pwm", data=data, background=bg(ppm(.Object)))
            return(pwm(.Object))
          })




setGeneric("ic",
           function(.Object) standardGeneric("ic"))
setMethod('ic',
          signature=signature(.Object="motif"),
          function(.Object){
            if(!is.null(.Object@ic))
              return(.Object@ic)
            if(is.null(ppm(.Object)))
              return(NULL)
            data <- data(ppm(.Object)) * data(pwm(.Object))
            data <- apply(data,2,sum)
            .ic(.Object)<-data
            return(ic(.Object))
          })


setGeneric("max.ic",
          function(.Object) standardGeneric("max.ic"))

setMethod("max.ic",
          signature=signature(.Object="motif"),
          function(.Object){
            return(log2(length(alphabet(.Object))))
          }
          )

        




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



