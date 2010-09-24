#!/usr/bin/Rscript

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

rd1 <- get(load(args[1]))
rd2 <-  get(load(args[2]))

space.rd1 <- levels(rd1$space)
space.rd2 <- levels(rd2$space)

#arg 1 is the query instance
#arg 2 is the subject, from which the nearest neighbours are found
res <- nearest(rd1,rd2)

#returns an integer vector containing the index of the nearest neighbour
#range in subject for each range in x.

ranges.rd1 <- ranges(rd1)
ranges.rd2 <-  ranges(rd2)

names(ranges.rd1) <- space.rd1
names(ranges.rd2) <- space.rd2

#if we don't have a matching space in the subject, we just don't look at that space.
res <- lapply(names(ranges.rd1), function(x){
     if(is.null(ranges.rd2[[x]])){
       warning(cat("No space '",x, "' in second (subject) ranged data object")) 
     }
     inds <- nearest(ranges.rd1[[x]], ranges.rd2[[x]])
     ranges.rd2[inds,]
   })

