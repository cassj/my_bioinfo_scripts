#!/usr/bin/Rscript

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

nearest <- get(load(args[1]))
file = args[2]
# file<-"xdnvev_v_restchip_tssplot.ps"

fc <- nearest[,"nearest.values.log2FoldChange"]
pval <- nearest[,"nearest.values.pVal"]
dist.tss <- nearest[,"dist.to.bindingsite"]

pval.cut <- 0.0001
do.plot <- function(file, xrange=NA ){
  postscript(file=file, paper="special", width=6, height=6)

  if(is.na(xrange))xrange <- max(abs(dist.tss))
  inds <- which(abs(dist.tss)<=xrange)
  plot(dist.tss[inds], fc[inds], pch=".")
  cords <- range(dist.tss[inds])

  inds <- intersect(which(pval <= pval.cut), inds)
  points(dist.tss[inds], fc[inds], col="red")

  lines(cords, c(1,1), col="blue")
  lines(cords, c(-1,-1), col="blue")
  dev.off()
}

files <-c(file,
          sub(".ps","10000kb.ps",file),
          sub(".ps","1000kb.ps",file),
          sub(".ps","100kb.ps",file),
          sub(".ps","10kb.ps",file),
          sub(".ps","5kb.ps",file),
          sub(".ps","1kb.ps",file))

xranges <- c(NA, 10000000, 1000000,100000,10000,5000, 1000)
for(i in 1:length(files)){
  do.plot(file=files[i], xrange=xranges[i])
}


#what are the genes with low p and low tss dist?
save.genes<-function(file, xrange=NA){
  file <- sub(".ps",".csv",file)
  
  if(is.na(xrange))xrange <- max(abs(dist.tss))
  inds <- which(abs(dist.tss) <= xrange)
  inds <- intersect(which(pval <= pval.cut), inds)

  my.cols <- c("dist.to.bindingsite",

               "start", "nearest.end", "values.neg10log10pVal",
               "values.feature", "values.distancetoFeature",
               "values.mgi_symbol", "values.description",

               "nearest.start", "nearest.end",
               "nearest.values.log2FoldChange", "nearest.values.pVal","nearest.values.FDR",
               "nearest.values.start_position", "nearest.values.end_position", "nearest.values.strand",
               "nearest.values.Symbols", "nearest.values.feature","nearest.values.mgi_symbol"
               )

  dat <- nearest[inds,my.cols]
  write.csv(dat, file=file)
  
}

for(i in 1:length(files)){
  save.genes(file=files[i], xrange=xranges[i])
}

