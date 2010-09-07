
qw <- function(...) {
  as.character(sys.call()[-1])
}


liftOver<-function(data, chain.file, ucsc.format=T, chr.col="chr", start.col="start",end.col="end"){

  #data should be a matrix or dataframe with cols for chr, start and end
  #TODO: Or a RangedData / IRange object, 
 
  this<-data.frame(chr=as.character(data[,chr.col]),start=as.numeric(data[,start.col]),end=as.numeric(data[,end.col]), stringsAsFactors=F)

  #Normal counting specifies ranges in 1-based, fully closed form.
  #UCSC specifies ranges in  0-based, half open
  if (!ucsc.format){
      this$start<-this$start-1
  }

  #use the chrstartend as an id so we can map back to the original data
  ids = paste(as.character(data[,chr.col]),data[,start.col],data[,end.col], sep=".")
  this <- cbind(this, ids)
  #If we have duplicate positions, remove them for the liftOver 
  this <- unique(this)
  #we also need to chuck out anything that doesn't have positional data.
  no.chr <- which(is.na(this[,1]))
  no.start <- which(is.na(this[,2]))
  no.end <- which(is.na(this[3]))

  inds <- unique(c(no.chr,no.start,no.end))

  if( length(inds)>0 ){ this <- this[-1*inds,] }
  
  
  ##all this stuff should be a .C() call but I don't have time to make it work just now.
  in.bed <- tempfile()

  #need to watch we don't get scientific notation printed out
  options(scipen=10)
  write.table(this, file=in.bed, sep="\t", row.names=F, col.names=F, quote=F)

  out.bed <- tempfile()
  out.um.bed <- tempfile()
  lo.cmd <- paste("liftOver", in.bed, chain.file, out.bed, out.um.bed)
  system(lo.cmd)

  try(
    new.bed<-read.table(out.bed, sep="\t") 
      ,silent=T)
  try(
    new.um <- read.table(out.um.bed, sep="\t")
      ,silent=T)
  
  #throw away the files
  unlink(c(in.bed, out.bed, out.um.bed))

  if (!exists('new.bed')){stop("No successful mappings")}

  #use the ids as rownames chuck out the column
  #order in the same order as our original data
  #which should stick NAs in for anything that can't be mapped
  rownames(new.bed) <- new.bed[,4]
  new.bed <- new.bed[,-4]
  new.bed <- new.bed[ids,]

  if(!ucsc.format){
   #put the data back to 1-based
   new.bed[,2] <- new.bed[,2]+1
  }

  #replace the new positions in the original dataset
  data[,chr.col] <- new.bed[,1]
  data[,start.col] <- new.bed[,2]
  data[,end.col] <- new.bed[,3]

  #TODO: return some information about the data that won't map
  
  return(data)
  
}



#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))


#function to do the liftover
liftOver.sorted <- function(data){
 
  #liftover and append to outfile
  #liftOver expects bed, which is chr start end and
  #intervals must be at least a base long
  to.map<-data.frame("chr"   = as.character(data[,"match.chr"]),
                     "start" = as.numeric(data[,"match.pos"])-1,
                     "end"   = as.numeric(data[,"match.pos"])
                     )

  mapped.mm9 <- liftOver(to.map, chain.file="lib/mm8ToMm9.over.chain")
  
  #and stick the results back in the dataframe
  data[,"match.chr"] <- as.character(mapped.mm9[,"chr"])
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

  #ditch anything where we don't have a chr
  ids <- which(is.na(new.data$match.chr))
  if (length(ids)>0){new.data <- new.data[(-1*ids),];}
  write.table(new.data, file=newfile, append=T, quote=F, sep="\t", row.names=F, col.names=F, na="")
}

close(in.file)


