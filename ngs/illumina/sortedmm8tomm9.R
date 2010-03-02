library(IRanges)
library(HadoopStreaming)
source('scripts/qw.R')
source('scripts/liftOver.R')

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))


#function to do the liftover
liftOver.sorted <- function(d, newfile ){

  #liftover and append to outfile
  #liftOver expects bed, which is chr start end and
  #intervals must be at least a base long
  to.map<-data[,qw(match.chr, match.pos, match.pos )]
  to.map[,2]<-to.map[,2]-1
                 
  colnames(to.map) <- c("chr","start","end")
  mapped.mm9 <- liftOver(to.map, chain.file="lib/mm8ToMm9.over.chain")
  
  #and stick the results back in the dataframe
  data[,qw(match.chr, match.pos)] <- mapped.mm9[,qw(chr,start)]
  data[,"match.pos"]<-data[,"match.pos"]+1

  #and append to the file
  write.table(data, file=newfile, append=T, quote=F, sep="\t", row.names=F, col.names=F)
}



#match pos +ve strand 1 based.
colnms<-qw(machine, run.number, lane, tile, x.coord, y.coord, index.string, read.number, read, quality.string,
         match.chr, match.contig, match.pos, match.strand, match.descriptor, single.read.aln.score,
         paired.read.aln.score, partner.chr, partner.contig, partner.offset, partner.strand)


#process the file as a stream, using the liftOver.sorted function
in.con <- textConnection(filename, open = "r")

newfile<-sub(".txt","_mm9.txt", filename)

#remove any existing newfile or we'll just keep appending to it.
unlink(newfile)

data <- hsTableReader(con,
                      col=colnms
                      chunkSize=100
                      FUN=liftOver.sorted(newfile=newfile),
                      ignoreKey=FALSE,
                      singleKey=FALSE)
close(con)



#for testing
data<-read.table(filename, nrow=100, sep="\t")
colnames(data)<-colnms

