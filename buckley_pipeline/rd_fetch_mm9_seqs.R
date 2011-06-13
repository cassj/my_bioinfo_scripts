#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ChIPpeakAnno)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]

#note that this gets the sequences for all of them, which
#could take a long time, depending.
rd <- get(load(filename))
rd.seqs = getAllPeakSequence(rd, upstream=0, downstream=0, genome=Mmusculus)

#save as a RangedData obj
new.filename <- sub('.R$', '_seqs.R', filename)
save(rd.seqs, file=new.filename)

#And as a csv file
new.filename <- sub('.R$', '_seqs.csv', filename)
df <- as.data.frame(rd.seqs)
write.csv(df, file=new.filename)


