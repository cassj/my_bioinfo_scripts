#!/usr/bin/Rscript

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

#for testing - why are the distances in bins?
#restchip = '../RESTChIP/Macs/NA_peaksWithSeqs.AnnoRangedData.R'
#xdnvev = '../XDNvEV/results/expression_data/AnnoRangedData_Limma.R'
 
query <- get(load(args[1]))
subject <-  get(load(args[2]))

#takes a query rd and a subject rd and returns
#a list with one entry per space containing a data.frame
#of the nearest ranges in subject to each of the ranges in
#query. Rows are named by query range name
nearest.rd <- function(query, subject ){
  spaces.query <- names(query)
  spaces.subject <- names(subject)
  common.spaces <- intersect(spaces.query, spaces.subject)
  if(length(common.spaces)<1){stop("No common spaces between subject and query")}

  res <- list()
  for(sp in common.spaces){
    cat("processing for space ",sp, "\n")
    q.iranges <- unlist(ranges(query[sp]))
    s.iranges <- unlist(ranges(subject[sp]))

   #nearest returns a integer vector containing the index of the 
   #nearest neightbour range in ‘subject’ for each range in ‘x’.
    inds <- nearest(q.iranges, s.iranges)
    
    #retrieve the info about nearest seqs
    nr <-  as.data.frame(subject[sp])[inds,]
    rownames(nr) <- sub(paste("^",sp,".",sep=""),'',names(q.iranges))
    colnames(nr) <- paste("nearest.", colnames(nr), sep="") 
    this.res <- cbind(as.data.frame(query[sp]),nr)

    #calculate the distance between the tss of the probe target and
    #the middle of the binding site.
    mid.bs <- apply(this.res[,c("nearest.start", "nearest.end")],1,mean)
    mid.bs <- round(mid.bs)
    genestart<-this.res[,"values.start_position"]
    inds <- which(this.res[,"values.strand"]==-1)
    genestart[inds] <- this.res[inds,"values.end_position"]
    dist.to.genestart<-genestart-mid.bs

    this.res <- cbind(dist.to.genestart,this.res)
    
    res[[sp]] <- this.res
  }
  return(res)  
}

res <- nearest.rd(query, subject)

#remove any spaces that didn't have any nearests 

#convert to a dataframe for all spaces
res.df <- do.call(rbind, res)
rownames(res.df) <- unlist(lapply(res, rownames))

filename = args[3]
save(res.df, file=paste(filename,".R",sep=""))

filename = paste(filename,".csv", sep="")
write.csv(res.df, file=filename)









# A note on the nearest function:
#     1. Find the ranges in ‘subject’ that overlap ‘xi’. If a single
#          range ‘si’ in ‘subject’ overlaps ‘xi’, ‘si’ is returned as
#          the nearest neighbor of ‘xi’. If there are multiple overlaps,
#          one of the overlapping ranges is chosen arbitrarily.
#
#       2. If no ranges in ‘subject’ overlap with ‘xi’, then the range
#          in ‘subject’ with the shortest distance from its end to the
#          start ‘xi’ or its start to the end of ‘xi’ is returned.
#
