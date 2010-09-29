#!/usr/bin/Rscript

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

nearest <- get(load(args[1]))
file = args[2]
#xrange <- as.numeric(args[4])

fc <- nearest[,"values.log2FoldChange"]
pval <- nearest[,"values.pVal"]
dist.tss <- nearest[,"dist.to.genestart"]

pval.cut <- 0.0001
do.plot <- function(nearest,file, xrange=NA ){
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
  do.plot(nearest=nearest, file=files[i], xrange=xranges[i])
}


#what are the genes with low p and low tss dist?
save.genes<-function(nearest, file, xrange=NA){
  file <- sub(".ps$",".csv",file)
  
  if(is.na(xrange))xrange <- max(abs(dist.tss))
  inds <- which(abs(dist.tss) <= xrange)
  inds <- intersect(which(pval <= pval.cut), inds)

  my.cols <- c("dist.to.genestart",
               "start", "end",
               "values.log2FoldChange", "values.pVal","values.FDR",
               "values.start_position", "values.end_position", "values.strand",
               "values.Symbols", "values.feature","values.mgi_symbol",
               "nearest.start", "nearest.end", "nearest.values.neg10log10pVal",
               "nearest.values.feature", "nearest.values.distancetoFeature",
               "nearest.values.mgi_symbol", "nearest.values.description")

  dat <- nearest[inds,my.cols]
  write.csv(dat, file=file)
  
}

for(i in 1:length(files)){
  save.genes(nearest=nearest, file=files[i], xrange=xranges[i])
}


#Multiple probes to a single gene. Is this giving misleading plot?
#take most \de probe?

files.nodup <-c(file,
                sub(".ps$","nodups.10000kb.ps",file),
                sub(".ps$","nodups.1000kb.ps",file),
                sub(".ps$","nodups.100kb.ps",file),
                sub(".ps$","nodups.10kb.ps",file),
                sub(".ps$","nodups.5kb.ps",file),
                sub(".ps$","nodups.1kb.ps",file))

unique.ensembl <- unique(nearest[,"values.feature"])
#sort by abs fc and just use the first instance of each gene.
sort.nr<-nearest[order(abs(nearest[,"values.log2FoldChange"]), decreasing=T),]

notdup<-!duplicated(sort.nr[,"values.feature"])
sort.nr<-sort.nr[notdup,]



for(i in 1:length(files)){
  do.plot(nearest=sort.nr, file=files.nodup[i], xrange=xranges[i])
  save.genes(nearest=sort.nr, file=files.nodup[i], xrange=xranges[i])

}

#doesn't make too much difference.
