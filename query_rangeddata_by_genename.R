#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("\nDESCRIPTION: Blah blah blah\n
\nUSAGE:\nquery_rd_by_genename.R <RangedDataFile> <genename> <genename> <genename>
OR:\nquery_rd_by_genename.R file <filename>\n\n ")
  options( show.error.messages=FALSE)
  stop()
}

library(IRanges)

#these files should only have 1 single R object in them
rdfile = args[1]
rd <- get(load(rdfile))

#and convert to dataframe as we're not doing overlap queries.
rd<-as.data.frame(rd)

if (args[2] == "file"){
  filename = args[3]
  if(length(filename) == 0){stop("No filename given?")}
  if(length(filename) > 1){stop("Can't handle multiple filenames")}
  genenames <- read.table(filename)[,1]
}else{
 genenames = args[-1]
}

#should probably check universe but not currently set.

inds <- lapply(genenames, function(x){rd[grep(x, rd[,"values.mgi_symbol"], ignore.case=T),]} )
names(inds) <- genenames

res <- rd[1,]

for(i in 1:length(inds)){
  this <- inds[[i]]
  if(nrow(this)>0){
    res<-rbind(res, this)
  }
}

res<-res[-1,]

write.csv(res)

