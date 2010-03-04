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
  to.map<-cbind("chr"   = data[,"match.chr"],
                "start" = as.numeric(data[,"match.pos"])-1,
                "end"   = as.numeric(data[,"match.pos"])
                )

  mapped.mm9 <- liftOver(to.map, chain.file="lib/mm8ToMm9.over.chain")
  
  #and stick the results back in the dataframe
  data[,qw(match.chr, match.pos)] <- mapped.mm9[,qw(chr,start)]
  data[,"match.pos"]<-data[,"match.pos"]+1

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
repeat {
  data <- t(as.data.frame(strsplit(readLines(in.file, n=1000, ok=T),"\t")))
  data<-cbind(data, rep("", nrow(data)))
  colnames(data) <- colnms
  if (nrow(data) == 0) break
  new.data <- liftOver.sorted(data)
  write.table(new.data, file=newfile, append=T, quote=F, sep="\t", row.names=F, col.names=F)
  
}

close(con)


#for testing
#data<-read.table(filename, nrow=100, sep="\t")
#colnames(data)<-colnms

