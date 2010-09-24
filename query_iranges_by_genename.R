#!/usr/local/bin/Rscript


args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("\nDESCRIPTION:\n
\nUSAGE:query_rangeddata_by_genename.R <RangedDataFile> <genename> <genename> <genename>
\nOR:query_rangeddata_by_genename.R file <RangedDataFile> <filename> ")
  options( show.error.messages=FALSE)
  stop()
}

args

#if (args[2] == "file"){
#  filename = args[3]
#  if(length(filename) == 0){stop("No filename given?")}
#  if(length(filename) > 1){stop("Can't handle multiple filenames")}
#  genenames <- read.data(
#}
#genenames = args
#
#genenames
