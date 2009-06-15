####
# This file contains the function make.ensembl.transcript.db()
# which takes an ensembl mart object and a filename as parameters.
# The mart must already have a dataset selected, so for example:
#
# ensmart <- useMart('ensembl_mart_46',
#                     dataset="mmusculus_gene_ensembl",
#                     archive=TRUE)
#
# make.ensembl.transcript.db(ensmart, 'my_database_file.db')
#
# The resulting sqlite database is created as 'my_database_file.db'
# the database file can be anything you like, but the format:
# species-genomebuild-featuretype.db is recommended
#
# The function returns TRUE on success
#
# Each chromosome gets a separate table, named Chr<whatever>
#
# Each table contains:
#
# Required fields for genomeLocation interface:
#   Chr - the chromosome name, in the form 'chrN'
#   GenomeStart - the start of the transcript in genome coordinates (+ve strand)
#   GenomeEnd - the end of the transcript in genome coordinates (+ve strand)
#   GenomeMid - the centre of the transcript in genome coords (+ve strand)
#   FeatureStart - the start of the transcript on the transcript strand
#   FeatureEnd - the end of the transcript on the transcript strand
#   Strand - the strand of the transcript
#   Description - Ensembl description of the transcript
# Extra Annotation fields:
#   EnsemblGeneID
#   EnsemblTranscriptID
#
###

library(biomaRt)
library(RSQLite)

make.ensembl.transcript.db <- function(ensmart=NULL, filename='transcripts.db'){

  #get the list of chromosome names in this mart.
  chrs <- getBM('chromosome_name', mart=ensmart)$"chromosome_name"
  if (length(chrs) < 1){stop("No chromosomes found for this dataset")}
  chrs <- chrs[-grep('NT_', chrs)]

  #get everything from the chromosome
  filters <- c('chromosome_name')
  attributes <- c('ensembl_transcript_id', 'description', 'ensembl_gene_id', 'transcript_start', 'transcript_end', 'strand', 'chromosome_name', 'transcript_status')
  
  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = filename)
  
  for(chr in chrs){
    
    cat("Generating for chr",chr,"\n");
    
    values<-list(chr)
    res <- getBM(attributes, filters=filters, values=values, mart=ensmart)
    f.start<-apply(res, 1,
               function(x){
                 if (x["strand"] == "-1"){
                   return(as.numeric(x["transcript_end"]))
               }else{
                 return(as.numeric(x["transcript_start"]))
               }
               }
               )
    f.end<-apply(res, 1,
               function(x){
                 if (x["strand"] == "-1"){
                   return(as.numeric(x["transcript_start"]))
               }else{
                 return(as.numeric(x["transcript_end"]))
               }
               }
               )
    #mean seems to have problems if the numbers are really big.
    mid <- as.numeric(res$'transcript_start') +
      
           round( ( as.numeric(res$'transcript_end') -
               as.numeric(res$'transcript_start')
              ) / 2
           )
    
    res<-data.frame(Chr=res$'chromosome_name',
                    GenomeStart=res$'transcript_start',
                    GenomeEnd=res$'transcript_end',
                    GenomeMid=mid,
                    FeatureStart=f.start,
                    FeatureEnd=f.end,
                    Strand=res$strand,
                    Description=res$description,
                    EnsemblGeneID=res$'ensembl_gene_id',
                    EnsemblTranscriptID=res$'ensembl_transcript_id')

    
    #sort by genome start
    res<-res[order(res$GenomeStart),]

    #make the table
    dbWriteTable(con, paste("chr",chr, sep=""), res)

    #create index on genome start, end, mid, feature start, end, 
    dbGetQuery(con,paste("CREATE INDEX 'genome_start_index_",chr,"' ON chr",chr,"('GenomeStart')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'genome_end_index_",chr,"' ON chr",chr,"('GenomeEnd')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'genome_mid_index_",chr,"' ON chr",chr,"('GenomeMid')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'feature_start_index_",chr,"' ON chr",chr,"('FeatureStart')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'feature_end_index_",chr,"' ON chr",chr,"('FeatureEnd')", sep=""));
    
    #and unique key on transcript ID
    dbGetQuery(con,paste("CREATE UNIQUE INDEX 'id_index_",chr,"' ON chr",chr,"('EnsemblTranscriptID')", sep=""))

  }

  dbDisconnect(con)
  return (TRUE);
}


make.ensembl.gene.db <- function(ensmart=NULL, filename='genes.db'){

  #get the list of chromosome names in this mart.
  chrs <- getBM('chromosome_name', mart=ensmart)$"chromosome_name"
  if (length(chrs) < 1){stop("No chromosomes found for this dataset")}
  chrs <- chrs[-grep('NT_', chrs)]

  #get everything from the chromosome
  filters <- c('chromosome_name')
  attributes <- c( 'ensembl_gene_id', 'description','start_position', 'end_position', 'strand', 'chromosome_name', 'status')
  
  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = filename)
  
  for(chr in chrs){
    
    cat("Generating for chr",chr,"\n");
    
    values<-list(chr)
    res <- getBM(attributes, filters=filters, values=values, mart=ensmart)
    f.start<-apply(res, 1,
               function(x){
                 if (x["strand"] == "-1"){
                   return(as.numeric(x["end_position"]))
               }else{
                 return(as.numeric(x["start_position"]))
               }
               }
               )
    f.end<-apply(res, 1,
               function(x){
                 if (x["strand"] == "-1"){
                   return(as.numeric(x["start_position"]))
               }else{
                 return(as.numeric(x["end_position"]))
               }
               }
               )
    #mean seems to have problems if the numbers are really big.
    mid <- as.numeric(res$'start_position') +
      
           round( ( as.numeric(res$'end_position') -
               as.numeric(res$'start_position')
              ) / 2
           )
    
    res<-data.frame(Chr=res$'chromosome_name',
                    GenomeStart=res$'start_position',
                    GenomeEnd=res$'end_position',
                    GenomeMid=mid,
                    FeatureStart=f.start,
                    FeatureEnd=f.end,
                    Strand=res$strand,
                    Description=res$description,
                    EnsemblGeneID=res$'ensembl_gene_id'
                    )

    
    #sort by genome start
    res<-res[order(res$GenomeStart),]

    #make the table
    dbWriteTable(con, paste("chr",chr, sep=""), res)

    #create index on genome start, end, mid, feature start, end, 
    dbGetQuery(con,paste("CREATE INDEX 'genome_start_index_",chr,"' ON chr",chr,"('GenomeStart')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'genome_end_index_",chr,"' ON chr",chr,"('GenomeEnd')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'genome_mid_index_",chr,"' ON chr",chr,"('GenomeMid')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'feature_start_index_",chr,"' ON chr",chr,"('FeatureStart')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'feature_end_index_",chr,"' ON chr",chr,"('FeatureEnd')", sep=""));
    
    #and unique key on gene ID
    dbGetQuery(con,paste("CREATE UNIQUE INDEX 'id_index_",chr,"' ON chr",chr,"('EnsemblGeneID')", sep=""))

  }

  dbDisconnect(con)
  return (TRUE);
}


make.ensembl.known.gene.db <- function(ensmart=NULL, filename='genes.db'){

  #get the list of chromosome names in this mart.
  chrs <- getBM('chromosome_name', mart=ensmart)$"chromosome_name"
  if (length(chrs) < 1){stop("No chromosomes found for this dataset")}
  chrs <- chrs[-grep('NT_', chrs)]

  #get everything from the chromosome
  filters <- c('chromosome_name','status')
  attributes <- c( 'ensembl_gene_id', 'description','start_position', 'end_position', 'strand', 'chromosome_name', 'status')
  
  m <- dbDriver("SQLite")
  con <- dbConnect(m, dbname = filename)
  
  for(chr in chrs){
    
    cat("Generating for chr",chr,"\n");
    
    values<-list(chr, 'known')
    res <- getBM(attributes, filters=filters, values=values, mart=ensmart)
    f.start<-apply(res, 1,
               function(x){
                 if (x["strand"] == "-1"){
                   return(as.numeric(x["end_position"]))
               }else{
                 return(as.numeric(x["start_position"]))
               }
               }
               )
    f.end<-apply(res, 1,
               function(x){
                 if (x["strand"] == "-1"){
                   return(as.numeric(x["start_position"]))
               }else{
                 return(as.numeric(x["end_position"]))
               }
               }
               )
    #mean seems to have problems if the numbers are really big.
    mid <- as.numeric(res$'start_position') +
      
           round( ( as.numeric(res$'end_position') -
               as.numeric(res$'start_position')
              ) / 2
           )
    
    res<-data.frame(Chr=res$'chromosome_name',
                    GenomeStart=res$'start_position',
                    GenomeEnd=res$'end_position',
                    GenomeMid=mid,
                    FeatureStart=f.start,
                    FeatureEnd=f.end,
                    Strand=res$strand,
                    Description=res$description,
                    EnsemblGeneID=res$'ensembl_gene_id'
                    )

    
    #sort by genome start
    res<-res[order(res$GenomeStart),]

    #make the table
    dbWriteTable(con, paste("chr",chr, sep=""), res)

    #create index on genome start, end, mid, feature start, end, 
    dbGetQuery(con,paste("CREATE INDEX 'genome_start_index_",chr,"' ON chr",chr,"('GenomeStart')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'genome_end_index_",chr,"' ON chr",chr,"('GenomeEnd')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'genome_mid_index_",chr,"' ON chr",chr,"('GenomeMid')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'feature_start_index_",chr,"' ON chr",chr,"('FeatureStart')", sep=""));
    dbGetQuery(con,paste("CREATE INDEX 'feature_end_index_",chr,"' ON chr",chr,"('FeatureEnd')", sep=""));
    
    #and unique key on transcript ID
    dbGetQuery(con,paste("CREATE UNIQUE INDEX 'id_index_",chr,"' ON chr",chr,"('EnsemblGeneID')", sep=""))

  }

  dbDisconnect(con)
  return (TRUE);
}
