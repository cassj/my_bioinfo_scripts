source('scripts/qw.R')
source('scripts/liftOver.R')

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))


#function to do the liftover
liftOver.sorted <- function(data){
 
  #liftover and append to outfile
  #liftOver expects bed, which is chr start end and
  #intervals must be at least a base long
  to.map<-data.frame("chr"   = data[,"match.chr"],
                     "start" = as.numeric(data[,"match.pos"])-1,
                     "end"   = as.numeric(data[,"match.pos"])
                     )

  mapped.mm9 <- liftOver(to.map, chain.file="lib/mm8ToMm9.over.chain")
  
  #and stick the results back in the dataframe
  data[,"match.chr"] <- mapped.mm9[,"chr"]
  data[,"match.pos"] <- 1+(mapped.mm9[,"start"] )

  return(data)
}


#Big file so process as a stream
in.file <- file(filename, "r")
newfile<-sub(".txt","_mm9.txt", filename)
unlink(newfile)

colnms<- qw(machine, run.number, lane, tile, x.coord, y.coord, index.string, read.number, read,
            quality.string, match.chr, match.contig, match.pos, match.strand, match.descriptor,
            single.read.aln.score, paired.read.aln.score, partner.chr, partner.contig,
            partner.offset, partner.strand)


## read in a chunk at a time, liftover and write to new file
n <- 1
repeat {
  cat(n,"\n")
  n <- n+1
  data <- NULL
  try({data <- read.table(in.file, sep="\t", header=F, nrow=10000, na.strings="")}, silent=T)
  if(is.null(data)) { break }
  colnames(data) <- colnms
  new.data <- liftOver.sorted(data)
  write.table(new.data, file=newfile, append=T, quote=F, sep="\t", row.names=F, col.names=F, na="")
}

close(in.file)


