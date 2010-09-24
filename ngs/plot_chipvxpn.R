#!/usr/bin/Rscript

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

#this should be nearest site to expression (so 1 row for each probe)
nearest <- get(load(args[1]))
file = args[2]
do.log <- args[3]
xrange <- args[4]

fc <- nearest[,"values.log2FoldChange"]
pval <- nearest[,"values.pVal"]
dist.tss <- nearest[,"nearest.values.distancetoFeature"]

if(!is.na(do.log) && do.log=="log") {
  neg <- which(dist.tss<0)
  dist.tss[neg] <- -1*dist.tss[neg] 
  dist.tss <- log(dist.tss)
  dist.tss[neg] <- -1*dist.tss[neg] 
}
cords <- range(dist.tss)
if(!is.na(xrange)){
  cords <- c((-1*xrange),xrange)
}

postscript(file=file, paper="special", width=6, height=6)
  plot(dist.tss, fc, pch=".", xlim=cords)
  inds <- which(pval<=0.001)
  points(dist.tss[inds], fc[inds], col="red", xlim=cords)
dev.off()
