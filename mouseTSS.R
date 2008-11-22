require(RSQLite)

# takes a vector of chromosome names (as chr1 etc),
# a vector (of the same length) of positions and
# a max number of transcripts to return for each position.
# Returns a list, containing $n and $transcripts eg:
#
# $transcripts
#         EnsemblGeneID EnsemblTranscriptID      TSS Distance  chr      pos
# 1  ENSMUSG00000051285  ENSMUST00000061280  7079231      258 chr1  7078973
# 2  ENSMUSG00000046101  ENSMUST00000052843  9898719    34392 chr1  9864327
# 3  ENSMUSG00000016918  ENSMUST00000088585 12708610   543352 chr1 12165258
# 4  ENSMUSG00000016918  ENSMUST00000088585 12708610    73744 chr1 12634866
# 5  ENSMUSG00000016918  ENSMUST00000088585 12708610    60999 chr1 12647611
# 6  ENSMUSG00000016918  ENSMUST00000088585 12708610    25110 chr1 12683500


#I wonder if it's worth doing this as a binary search?
#it would avoid having to calc TSS-pos for every entry in the
#table each time 

#Retrieve the closest Ensembl KNOWN gene
get.closest.tss <- function(chr,pos, n=1){
  if (length(chr)!=length(pos)){stop("chr must be the same length as pos")}

  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = "/home/cassj/work/MillHill/MASH1/MusMusculus_TSS.db")

  res <- data.frame()

  for(i in 1:length(chr)){
    q <- paste('select ', chr[i],'.*, abs("TSS"-',pos[i],') as Distance from ',chr[i],' order by Distance asc limit ',n, sep="")
    stuff <- dbGetQuery(con,q)
    stuff <- cbind(stuff[,-1], BS.Centre=pos[i])
    res<-rbind(res,stuff) 

  }

  dbDisconnect(con)
  return(list(transcripts=res,n=n))
}


#Retrieve the closest Ensembl KNOWN and Predicted genes
get.closest.predicted.tss <- function(chr,pos, n=1){
  if (length(chr)!=length(pos)){stop("chr must be the same length as pos")}

  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = "/home/cassj/work/MillHill/MASH1/MusMusculus_TSS_Pred.db")

  res <- data.frame()

  for(i in 1:length(chr)){
    q <- paste('select ', chr[i],'.*, abs("TSS"-',pos[i],') as Distance from ',chr[i],' order by Distance asc limit ',n, sep="")
    stuff <- dbGetQuery(con,q)
    stuff <- cbind(stuff[,-1], BS.Centre=pos[i])
    res<-rbind(res,stuff) 

  }

  dbDisconnect(con)
  return(list(transcripts=res,n=n))
}


#retrieve the closest RefSeq Peptide Gene
get.closest.refseq <- function(chr,pos, n=1){
  if (length(chr)!=length(pos)){stop("chr must be the same length as pos")}

  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = "/home/cassj/work/MillHill/MASH1/MusMusculus_TSS_RefSeq.db")

  res <- data.frame()

  for(i in 1:length(chr)){
    q <- paste('select ', chr[i],'.*, abs("TSS"-',pos[i],') as Distance from ',chr[i],' order by Distance asc limit ',n, sep="")
    stuff <- dbGetQuery(con,q)
    stuff <- cbind(stuff[,-1], BS.Centre=pos[i])
    res<-rbind(res,stuff) 

  }

  dbDisconnect(con)
  return(list(transcripts=res,n=n))
}





#Retrieve flanking Ensembl KNOWN gene
get.flanking.tss <- function(chr,pos, n=1){
  if (length(chr)!=length(pos)){stop("chr must be the same length as pos")}

  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = "/home/cassj/work/MillHill/MASH1/MusMusculus_TSS.db")

  res <- data.frame()

  for(i in 1:length(chr)){
    q <- paste('select ', chr[i],'.*, abs("TSS"-',pos[i],') as Distance from ',chr[i],' order by Distance asc limit ',n, sep="")
    stuff <- dbGetQuery(con,q)
    stuff <- cbind(stuff[,-1], BS.Centre=pos[i])
    res<-rbind(res,stuff) 

  }

  dbDisconnect(con)
  return(list(transcripts=res,n=n))
}



#retrive flanking RefSeq genes.
get.flanking.refseq <- function(chr,pos, n=1){
  if (length(chr)!=length(pos)){stop("chr must be the same length as pos")}

  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = "/home/cassj/work/MillHill/MASH1/MusMusculus_TSS_RefSeq.db")

  res <- data.frame()

  for(i in 1:length(chr)){
    q <- paste('select ', chr[i],'.*, abs("TSS"-',pos[i],') as Distance from ',chr[i],' order by Distance asc limit ',n, sep="")
    stuff <- dbGetQuery(con,q)
    stuff <- cbind(stuff[,-1], BS.Centre=pos[i])
    res<-rbind(res,stuff) 

  }

  dbDisconnect(con)
  return(list(transcripts=res,n=n))
}
