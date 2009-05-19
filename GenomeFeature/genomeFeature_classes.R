require(RSQLite)

### TODO:
#when creating databases of gFs, if a chr has been investigated and no
#features have been found on it, then it should get an empty table
#in the database
#
#Also in the get.nearest &c functions, we need to deal with the case of
#a chromosome not being in the database - ie not having been investigated,
#
# actually - this will all be resolved in postgres. 


#we're likely to have to deal with large lists of
#genomeFeatures, so they should be a dataframe
setClass("genomeFeature",
         representation(
                        features = "matrixy",
                        data.source="character",
                        species= "character",
                        genomeBuild="character",
                        creation.date="character",
                        notes="character",
                        feature.id.field="character",
                        name="character",
                        description="character"
                        ),

         )


setValidity("genomeFeature",
            function(object){
              nms<-colnames(object@features)
              if(!('Chr' %in% nms)){
                return("features must contain a Chr column")
              }
              if(!('GenomeStart' %in% nms)){
                return("features must contain a GenomeStart column")
              }
              if(!('GenomeEnd' %in% nms)){
                return("features must contain a GenomeEnd column")
              }
              if(!('Strand' %in% nms)){
                return("features must contain a Strand column")
              }
              if(!(object@feature.id.field == "")){
                if (!(object@feature.id.field %in% nms)){
                  return("feature.id.field not found in colnames of features")
                }
              }
              if(sum(duplicated(object@features[,object@feature.id.field]))!=0){
                return("Duplicated feature IDs are not allowed. If you haven't specified an ID field, this means you've got duplicate genome positions, either remove the duplicates or add a unique ID column")
              }
              return(TRUE);
            }
            )


#and set calculated cols if not already set
setMethod("initialize",
          "genomeFeature",
          function(.Object,
                   features=data.frame(
                     Chr = character(),
                     GenomeStart = numeric(),
                     GenomeEnd = numeric(),
                     Strand = numeric(),
                     Description=character(),
                     stringsAsFactors = FALSE
                     ),
                   data.source='',
                   species='',
                   genomeBuild='',
                   creation.date=date(),
                   notes='',
                   feature.id.field='',
                   name="",
                   description=""
                   )
          {
            
            nms<-colnames(features)
            
                                        #calculate other cols if they're not already there.
            if(!('GenomeMid' %in% nms)){
              GenomeMid <- apply(cbind(features[,"GenomeStart"],
                                       features[,"GenomeEnd"]), 1,mean)
              GenomeMid <- round(GenomeMid)
              features <- cbind(features, GenomeMid=GenomeMid)
            }
            
            if(!('FeatureStart' %in% nms)){
              FeatureStart <- features[,"GenomeStart"]
              inds <- which(features[,"Strand"] == -1)
              FeatureStart[inds] <- features[inds,"GenomeEnd"]
              features <- cbind(features, FeatureStart=FeatureStart)
              
            }
            if(!('FeatureEnd' %in% nms)){
              FeatureEnd <- features[,"GenomeEnd"]
              inds <- which(features[,"Strand"]==-1)
              FeatureEnd[inds] <- features[inds,"GenomeStart"]
              features <- cbind(features, FeatureEnd=FeatureEnd)
            }
            if(feature.id.field==""){
              feature.id.field="ID"
              ID<-paste(features[,"Chr"], paste(features[,"GenomeStart"], features[,"GenomeEnd"], sep="-"), sep=":")
              features <- cbind(ID, features)
            }

            #force chr to be char,  not factor
            features[,"Chr"]<-as.character(features[,"Chr"])            
            features[,feature.id.field]<- as.character(features[,feature.id.field])
            
            .Object@features <- features
            .Object@data.source <- data.source
            .Object@species <- species
            .Object@genomeBuild <- genomeBuild
            .Object@creation.date <- creation.date
            .Object@notes <- notes;
            .Object@feature.id.field <- feature.id.field
            .Object@name <- name
            .Object@description <- description
            
            validObject(.Object)
            return(.Object)
          }
          )


# user-friendly constructors

#basic constructor, just passes all args to new
gF <- genomeFeature <- function(...){
  new ("genomeFeature", ...)
}

#constructor from bed file
#name arg overrides header name
gF.from.bed <-genomeFeature.from.bed <- function(filename=NULL,
                                                 head=TRUE,
                                                 sep=" ",
                                                 description="",
                                                 notes="",
                                                 name = "",
                                                 ...){
  if(is.null(filename)){
    stop("No filename given")
  }
  else{
    skip=0
    if(head){
      header <- readLines(filename, n=1)
      if(notes == ""){
        notes <- paste("BED Header:", header, sep=" ")
        if (name=="") {name <- gsub('.*name=[\',"]?(\\w+)\\S.*','\\1',header, perl=TRUE)}
        if(description==""){description <- gsub('.*description=[\',"](.*)[\',"].*','\\1',header, perl=TRUE)}
      }
      skip=1
    }
    bed<-read.csv(file=filename, skip=skip, header=FALSE, sep=sep)
    cols <- c("Chr","GenomeStart","GenomeEnd", "Score", "Strand")
    colnames(bed) <- cols[1:ncol(bed)]
    if(!('Strand' %in% colnames(bed))){
      bed <- cbind(bed, Strand=rep(0, nrow(bed)))
    }
  }
  
  new ("genomeFeature", features=bed, data.source=filename, notes=notes, name=name, description=description, ...)
}


#constructor from database
gF.from.db <-genomeFeature.from.db <- function(filename=NULL,
                                               data.source="",
                                               notes="",
                                               name="",
                                               description="",
                                               ...){
  if(is.null(filename)){
    stop("No filename given")
  }
  else{
    if(name != ""){name=filename}
    m <- dbDriver("SQLite")
    con <- dbConnect(m, dbname = filename)
    tables <- dbListTables(con)
    data=NULL
    for (table in tables){
      these <-  dbGetQuery(con,paste("SELECT * FROM ",table, sep=""));
      if(!is.null(data)){ data <- rbind(data, these)}
      else{data <- these}
    } 
  }
  
  new ("genomeFeature", features=data, data.source=filename, notes=notes, name=name, description=description, ...)
}











# gets


if (!isGeneric("features")) {
     if (is.function("features"))
       fun <- features
     else fun <- function(object) standardGeneric("features")
     setGeneric("features", fun)
}

setMethod("features", "genomeFeature", function(object) object@features)


if (!isGeneric("data.source")) {
     if (is.function("data.source"))
       fun <- data.source
     else fun <- function(object) standardGeneric("data.source")
     setGeneric("data.source", fun)
}

setMethod("data.source", "genomeFeature", function(object) object@data.source)


if (!isGeneric("species")) {
     if (is.function("species"))
       fun <- species
     else fun <- function(object) standardGeneric("species")
     setGeneric("species", fun)
}

setMethod("species", "genomeFeature", function(object) object@species)


if (!isGeneric("genomeBuild")) {
     if (is.function("genomeBuild"))
       fun <- genomeBuild
     else fun <- function(object) standardGeneric("genomeBuild")
     setGeneric("genomeBuild", fun)
}

setMethod("genomeBuild", "genomeFeature", function(object) object@genomeBuild)

if (!isGeneric("creation.date")) {
     if (is.function("creation.date"))
       fun <- creation.date
     else fun <- function(object) standardGeneric("creation.date")
     setGeneric("creation.date", fun)
}

setMethod("creation.date", "genomeFeature", function(object) object@creation.date)

if (!isGeneric("notes")) {
     if (is.function("notes"))
       fun <- notes
     else fun <- function(object) standardGeneric("notes")
     setGeneric("notes", fun)
}

setMethod("notes", "genomeFeature", function(object) object@notes)

if (!isGeneric("feature.id.field")) {
     if (is.function("feature.id.field"))
       fun <- feature.id.field
     else fun <- function(object) standardGeneric("feature.id.field")
     setGeneric("feature.id.field", fun)
}

setMethod("feature.id.field", "genomeFeature", function(object) object@feature.id.field)



# sets

#Once the features are set, you can't change them, on the grounds that
#if you're changing the features, you should really be updating every
#other slot appropriately. So create a new instance.
#setGeneric("features<-", function(x, value) standardGeneric("features<-"))
#setReplaceMethod("features", "genomeFeature", function(x, value) {
#  x@features <- value
#  x
#})


setGeneric("data.source<-", function(x, value) standardGeneric("data.source<-"))
setReplaceMethod("data.source", "genomeFeature", function(x, value) {
  x@data.source <- value
  x
})


setGeneric("species<-", function(x, value) standardGeneric("species<-"))
setReplaceMethod("species", "genomeFeature", function(x, value) {
  x@species <- value
  x
})


setGeneric("genomeBuild<-", function(x, value) standardGeneric("genomeBuild<-"))
setReplaceMethod("genomeBuild", "genomeFeature", function(x, value) {
  x@genomeBuild <- value
  x
})


setGeneric("creation.date<-", function(x, value) standardGeneric("creation.date<-"))
setReplaceMethod("creation.date", "genomeFeature", function(x, value) {
  x@creation.date <- value
  x
})


setGeneric("notes<-", function(x, value) standardGeneric("notes<-"))
setReplaceMethod("notes", "genomeFeature", function(x, value) {
  x@notes <- value
  x
})


setGeneric("feature.id.field<-", function(x, value) standardGeneric("feature.id.field<-"))
setReplaceMethod("feature.id.field", "genomeFeature", function(x, value) {
  x@notes <- value
  x
})





#method to write as db

if (!isGeneric("make.feature.db")) {
     if (is.function("make.feature.db"))
       fun <- features
     else fun <- function(object, filename) standardGeneric("make.feature.db")
     setGeneric("make.feature.db", fun)
}

#Note to self: need to fix this to take into account the possibility of missing chrms

setMethod("make.feature.db", "genomeFeature",
          
   function(object, filename="features.db"){

     m <- dbDriver("SQLite")
     con <- dbConnect(m, dbname = filename)
     chrs <- unique(features(object)[,"Chr"])
     feats <- features(object)
     
     for(chr in as.character(chrs)){
                                        #get just this chr.
       res <- feats[feats[,"Chr"]==chr,]
                                        #sort by GenomeEnd
       res<-res[order(res[,"GenomeEnd"]),]
       
                                        #sort by GenomeStart
       res<-res[order(res[,"GenomeStart"]),]
       
                                        #make the table
       dbWriteTable(con, chr, res)
       
                                        #create index on genome start, end, mid, feature start, end, 
       dbGetQuery(con,paste("CREATE INDEX 'genome_start_index_",chr,"' ON ",chr,"('GenomeStart')", sep=""));
       dbGetQuery(con,paste("CREATE INDEX 'genome_end_index_",chr,"' ON ",chr,"('GenomeEnd')", sep=""));
       dbGetQuery(con,paste("CREATE INDEX 'genome_mid_index_",chr,"' ON ",chr,"('GenomeMid')", sep=""));
       dbGetQuery(con,paste("CREATE INDEX 'feature_start_index_",chr,"' ON ",chr,"('FeatureStart')", sep=""));
       dbGetQuery(con,paste("CREATE INDEX 'feature_end_index_",chr,"' ON ",chr,"('FeatureEnd')", sep=""));

       
                                        #and if we have a unique ID field, create a unique index on that
       if(feature.id.field(object)!=''){
         q<-paste("CREATE UNIQUE INDEX 'id_index_",chr,"' ON ",chr,"('",feature.id.field(object),"')", sep="")
         dbGetQuery(con,q)
       }
     }
     dbDisconnect(con)
     return(TRUE)
   }
          )







#find nearest n features 

if (!isGeneric("get.nearest")) {
     if (is.function("get.nearest"))
       fun <- get.nearest
     else fun <- function(object, compare.to, field.1="GenomeMid", field.2="GenomeMid",n=1, feature.context=0, just.data.2=FALSE, same.strand=TRUE ) standardGeneric("get.nearest")
     setGeneric("get.nearest", fun)
}


##if compare.to is a genomeFeature, make a tmp db and use that
#will implement in just R later.
setMethod("get.nearest", c("genomeFeature", "genomeFeature"),
          function(object, compare.to, field.1="GenomeMid", field.2="GenomeMid", n=1, feature.context=0, just.data.2=FALSE, same.strand=TRUE){
            filename=tempfile()
            make.feature.db(compare.to, filename=filename)
            get.nearest(object, compare.to=filename,field.1=field.1, field.2=field.2, n=n, feature.context=feature.context, just.data.2=just.data.2 )
            
          })

       
#and if compere.to is character, interpret as db filename
setMethod("get.nearest", c("genomeFeature","character"),
          
          function(object, compare.to, field.1="GenomeMid", field.2="GenomeMid", n=1, feature.context=0, just.data.2=FALSE, same.strand=TRUE){
            
                                        #note to self- this would be much better if we used SQlite's
                                        #spatialite thing.
            m <- dbDriver("SQLite")
            con <- dbConnect(m, dbname = compare.to)
            
            feats <- features(object)
            chr <- as.character(feats[,"Chr"])
            pos <- as.character(feats[,field.1])
            res <- data.frame()
            colnames(feats) <- paste("feat.", colnames(feats), sep="")
            strand <- feats[,"feat.Strand"]
            
            
            for(i in 1:nrow(feats)){
              cat(i,"\n");
              
              nearest<-NA
              these <- feats[i,]
              
              if (dbExistsTable(con, chr[i])){
                if (same.strand & feats[i,"feat.Strand"]!=0){
                  q <- paste('select ', chr[i],'.*, abs("',field.2,'"-',pos[i],') as Distance from ',chr[i],' where Strand=0 or Strand=',strand[i],' order by Distance asc limit ',n, sep="")
                }else{
                  q <- paste('select ', chr[i],'.*, abs("',field.2,'"-',pos[i],') as Distance from ',chr[i],' order by Distance asc limit ',n, sep="")
                }
                nearest <- dbGetQuery(con,q)
              }else{

                q <- paste('select chr1.*, abs(\"FeatureStart\") as Distance from chr1 limit 1', sep="")
                nearest <- dbGetQuery(con,q)
                nearest[1:length(nearest)] <- NA
              }
              
              nearest <- cbind(nearest[,-1]) #get rid of the row.names
              if(nrow(nearest)==0){next()}
              colnames(nearest)<-paste("nearest.",colnames(nearest), sep="")
              
              data<-parse.nearest(these,
                                  nearest,
                                  paste("feat.", field.1, sep=""),
                                  paste("nearest.", field.2, sep=""),
                                  feature.context,
                                  just.data.2,
                                  same.strand)
              
              
              res<-rbind(res, data )
              
            }
            
            dbDisconnect(con)
            
                                        #create a comparison obj
      res<-nF(data=res,
              species=species(object),
              genomeFeature="todo",
              feature.context=feature.context,
              n=n,
              just.data.2=just.data.2,
              field.1=field.1,
              field.2=field.2)
            return(res)
          })



parse.nearest <- function(these, nearest, field.1, field.2, feature.context, just.data.2, same.strand){

  
  
    if(nrow(nearest)>1){
      for(i in 2:nrow(nearest)){
        these <- rbind(these,these[1,])
      }
    }
    
    
    if (feature.context == 0){
      data <- nearest
    }
    
                                        #if feature.context=1, context is feature & feature strand
                                        #if feature strand is 0, use top strand (NOT nearest strand)
    if(feature.context==1){
      down.ids<-which(nearest[,field.2] <= these[,field.1])
      up.ids<-which(nearest[,field.2] > these[,field.1])
      
      orientation<-rep('UP',nrow(nearest))
      orientation[down.ids]<-"DOWN"
                                        #and invert -ve strand features 
      neg<-which(these[,"feat.Strand"]==-1)
      orientation[intersect(down.ids,neg)]="UP"
      orientation[intersect(up.ids, neg)]="DOWN"
    
      opp <- intersect(nearest[,"nearest.Strand"]!=0, intersect(these[,"feat.Strand"]!=0, these[,"feat.Strand"]!=nearest[,"nearest.Strand"]))
      orientation[opp] <- paste(orientation[opp], "OPP", sep=".")
      data <- cbind(nearest, orientation)
      
    }
                                        #if feature.context=2, context is nearest & nearest strand
                                        #if nearest strand is 0, use top strand (NOT feature strand)
    if(feature.context==2){
      up.ids<-which(nearest[,field.2] > these[,field.1])
      down.ids<-which(nearest[,field.2] <= these[,field.1])
      
      orientation<-rep('UP',nrow(nearest))
      orientation[down.ids]<-"DOWN"
                                        #and invert -ve strand features 
      neg<-which(nearest[,"nearest.Strand"]==-1)
      orientation[intersect(down.ids,neg)]="UP"
    orientation[intersect(up.ids, neg)]="DOWN"
      
    
      opp <- intersect(these[,"feat.Strand"]!=0, intersect(nearest[,"nearest.Strand"]!=0, these[,"feat.Strand"]!=nearest[,"nearest.Strand"]))
      orientation[opp] <- paste(orientation[opp], "OPP", sep=".")
      
      data <- cbind(nearest, orientation)
    }
  
  if(!just.data.2){ data <- cbind(these, data)}
  
}






#basically the same as get.nearest, returns an object of class featureNearest, but
#returns n features up- and n features down- stream

# fix this to work with new parse.nearest
if (!isGeneric("get.flanking")) {
     if (is.function("get.flanking"))
       fun <- get.flanking
     else fun <- function(object, compare.to, field.1="GenomeMid", field.2="GenomeMid", feature.context=FALSE, just.data.2=FALSE, n=1 ) standardGeneric("get.flanking")
     setGeneric("get.flanking", fun)
}




##if compare.to is a genomeFeature, make a tmp db and use that
#will implement in just R later.
setMethod("get.flanking", c("genomeFeature", "genomeFeature"),
          function(object, compare.to, field.1="GenomeMid", field.2="GenomeMid", feature.context=FALSE, just.data.2=FALSE, n=1){
            filename=tempfile()
            make.feature.db(compare.to, filename=filename)
            get.flanking(object, compare.to=filename,field.1=field.1, field.2=field.2, n=n, feature.context=feature.context, just.data.2=just.data.2 )
            
          })

setMethod("get.flanking", c("genomeFeature","character"),
          function(object, compare.to, field.1="GenomeMid", field.2="GenomeMid", feature.context=FALSE, just.data.2=FALSE, n=1 ){

            m <- dbDriver("SQLite")
            con <- dbConnect(m, dbname = compare.to)

            feats <- features(object)
            chr <- as.character(feats[,"Chr"])
            pos <- as.character(feats[,field.1])
            res <- data.frame()
            strand <- as.numeric(feats[,"Strand"])
            colnames(feats) <- paste("feat.", colnames(feats), sep="")

            for(i in 1:nrow(feats)){
              
              qUp <- paste('select ', chr[i],'.*, "',field.2,'"-',pos[i],' as Distance from ',chr[i],' where Distance > 0 order by Distance asc limit ',n, sep="")
              qDown <- paste('select ', chr[i],'.*, "',pos[i],'"-',field.2,' as Distance from ',chr[i],' where Distance>0 order by Distance asc limit ',n, sep="")

              up <- dbGetQuery(con, qUp)
              down <- dbGetQuery(con, qDown)
              up <- up[,-1]
              down <- down[,-1]
              nearest <- rbind(up, down)
              if(nrow(nearest)==0){next()}
              colnames(nearest)<-paste("nearest.",colnames(nearest), sep="")
              these <- feats[i,]


              data<-parse.nearest(these,
                                  nearest,
                                  paste("feat.", field.1, sep=""),
                                  paste("nearest.", field.2, sep=""),
                                  feature.context,
                                  just.data.2)
   

              res<-rbind(res, data )
              
            }
            
            
            dbDisconnect(con)
            
                                        #create a comparison obj
            res<-nF(data=res,
                    species=species(object),
                    genomeFeature="todo",
                    feature.context=feature.context,
                    n=n,
                    just.data.2=just.data.2,
                    field.1=field.1,
                    field.2=field.2)
            return(res)
          })









#Output as a bed file

#Need to have a think about how we deal with IDs. 
if (!isGeneric("make.bed")) {
     if (is.function("make.bed"))
       fun <- make.bed
     else fun <- function(object, filename="features.bed", header=TRUE, track.name="Features", track.description="A list of features", track.score.field=NULL ) standardGeneric("make.bed")
     setGeneric("make.bed", fun)
}

setMethod("make.bed", "genomeFeature",
          
   function(object, filename="features.bed", header=TRUE, track.name="Features",track.description="A list of features", track.score.field=NULL){
     
     if(header){
       cat(c("track",paste("name='",track.name, "'",sep=""), paste("description='", track.description,"'", sep="")), file=filename, sep="\t")
       if(!is.null(track.score.field)){
         cat("\tscore=1", file=filename,append=TRUE)
        }
       cat("\n",file=filename, append=TRUE)

     }
     
     if(!is.null(track.score.field)){
       mat <- object@features[,c("Chr", "GenomeStart","GenomeEnd", track.score.field)]

     }else{
       mat <- object@features[,c("Chr", "GenomeStart","GenomeEnd")]
     }

     #make the ID chr_start_end
    if (!is.null(object@feature.id.field)){
      feat.ids <-  object@features[,c(object@feature.id.field)]
    }else{
      feat.ids <- paste(mat[,1],mat[,2],mat[,3], sep="_")
    }
      
     mat <- cbind(mat,feat.ids)
     write.table(mat, file=filename, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
     
   }
 )












# get within: retrieve any features which are within up bases upstream and
# down bases downstream of the features in object. Setting up and down =0
# (the default) gets overlapping features.
# up and down are by default interpreted in the genome (+ve strand) context
# (feature.context=0) using GenomeStart and GenomeEnd positions for ref
# setting feature.context to 1 will mean they are interpreted in the context of
# the feature strand, using FeatureStart and FeatureEnd as a reference
# Unlike get.nearest, feature.context=2 is meaningless.
# if same.strand=TRUE, then only features on the same strand will be retrieved.
# If feature strand is 0 then genome context will be used regardless of the
# settings for feature.context and same.strand will be ignored
# similarly, if same.strand=TRUE, found features with strand=0 will be accepted

if (!isGeneric("get.within")) {
     if (is.function("get.within"))
       fun <- get.within
     else fun <- function(object, compare.to, feature.context=0, just.data.2=FALSE, up=0, down=0, same.strand=TRUE ) standardGeneric("get.within")
     setGeneric("get.within", fun)
}



#for now, if we get a genomeFeature, make a DB of it and use that
setMethod("get.within", c("genomeFeature", "genomeFeature"),
          function(object, compare.to, feature.context=0, just.data.2=FALSE, up=0, down=0, same.strand=TRUE){
            filename=tempfile()
            make.feature.db(compare.to, filename=filename)
            get.within(object, compare.to=filename, feature.context=feature.context, just.data.2=just.data.2, up=up, down=down, same.strand=same.strand )
            
          })

setMethod("get.within", c("genomeFeature","character"),
          function(object, compare.to, feature.context=0, just.data.2=FALSE, up=0, down=0, same.strand=TRUE ){

            m <- dbDriver("SQLite")
            con <- dbConnect(m, dbname = compare.to)

            feats <- features(object)
            chr <- as.character(feats[,"Chr"])
            res <- data.frame()
            strand <- as.numeric(feats[,"Strand"])
            start<-feats[,'GenomeStart']-up
            end <- feats[,'GenomeEnd']+down
            res <- data.frame()
            colnames(feats) <- paste("feat.", colnames(feats), sep="")

            for(i in 1:nrow(feats)){
              cat(i, 'of', nrow(feats), "\n")
              if(same.strand & feats[i,"feat.Strand"]!=0){
                q1<-paste('select ', chr[i],'.* from ', chr[i], ' where GenomeStart >',start[i],' and GenomeEnd < ',end[i],, 'and (Strand=0 or Strand=',strand[i], ')', sep="")
                q2<-paste('select ', chr[i],'.* from ', chr[i], ' where GenomeEnd >',start[i],' and GenomeStart < ',end[i],,' and (Strand=0 or Strand=',strand[i], ')',  sep="")
                
              }else{
                q1<-paste('select ', chr[i],'.* from ', chr[i], ' where GenomeStart >',start[i],' and GenomeEnd < ',end[i], sep="")
                q2<-paste('select ', chr[i],'.* from ', chr[i], ' where GenomeEnd >',start[i],' and GenomeStart < ',end[i], sep="")
              }

              ol1 <- dbGetQuery(con, q1)
              ol2 <- dbGetQuery(con, q2)
              ol1 <- ol1[,-1]
              ol2 <- ol2[,-1]
              overlapping <- rbind(ol1, ol2)
              if(nrow(overlapping)==0){next()}
              colnames(overlapping)<-paste("overlapping.",colnames(overlapping), sep="")
              these <- feats[i,]
               
              if(nrow(overlapping)>1){
                for(i in 2:nrow(overlapping)){
                  these <- rbind(these,these[1,])
                }
              }
              these <- cbind(these,overlapping)
              res<-rbind(res, these)
            }

            
            dbDisconnect(con)
            
                                        #create a comparison obj
            res<-oF(data=res,
                    species=species(object),
                    genomeFeature="todo",
                    feature.context=feature.context,
                    just.data.2=just.data.2,
                    up=up,
                    down=down
                    )
            return(res)
          })








#### just some useful functions (not methods) ###

# you need to use biomaRt to create an ensmart

# ensmart <- useMart('ensembl_mart_46',
#                     dataset="mmusculus_gene_ensembl",
#                     archive=TRUE)
#
# and nearest should be the result of calling the get.nearest
# method with an ensembl_transcript database.

get.ensembl.transcript.context <- function(nearest, ensmart){


  #change this when you've got some accessor methods
  nearest <- nearest@data
  exon <- utr <- intron <- as.character(rep(NA,nrow(nearest)))

  #only need context for intragenic features:
  intra.ids <- which( ( nearest$feat.GenomeMid > nearest$nearest.GenomeStart )
                      & ( nearest$feat.GenomeMid<nearest$nearest.GenomeEnd ) )

  
  #ok, grab those transcripts from ensembl.
  filters <- c('ensembl_gene_id')
  values <- nearest[intra.ids,"nearest.EnsemblGeneID"]
  attributes <- c(
                  'sequence_exon_stable_id',
                  'structure_exon_chrom_start',
                  'structure_exon_chrom_end',
                  'structure_exon_rank',
                  'sequence_exon_chrom_strand',
                  'sequence_3utr_start',
                  'sequence_3utr_end',
                  'sequence_5utr_start',
                  'sequence_5utr_end',
                  'ensembl_transcript_id',
                  'ensembl_gene_id'
                 )

  
#  #above doesn't work for mm8.
#  attributes<-c('sequence_exon_stable_id',
#                'sequence_exon_chrom_start',
#                'sequence_exon_chrom_end',
#                'rank',
#                'sequence_exon_chrom_strand',
#                '3_utr_start',
#                '3_utr_end',
#                '5_utr_start',
#                '5_utr_end',
#                'ensembl_transcript_id',
#                'ensembl_gene_id'
#                )

  res <- getBM(attributes, filters=filters, values=values, mart=ensmart)

  n<-0
  len<-length(intra.ids)
  #ok, now for each transcript ID in intra.ids, check where the genome pos is
  for(i in intra.ids){
    cat(i,"\n")
    n<-n+1
    cat(n,'of',len,"\n" )
    ts<-nearest$nearest.EnsemblGeneID[i]
    these <- res[res$'ensembl_gene_id'== ts,]
    pos <- nearest$feat.GenomeMid[i]
    start.ok <- which(these$'exon_chrom_start'<=pos)
    end.ok <- which(these$'exon_chrom_end'>= pos)
    ex <- intersect(start.ok, end.ok)
    
    if(length(ex)!=0){
                                           #exonic,
      exon[i] <- paste(these[ex,'sequence_exon_stable_id'],' (' ,these[ex,'structure_exon_rank'], ')',sep="")
                                        #see if it's in the utr
      if(!is.na(these[ex, 'sequence_3utr_start'])){
        if(pos>=these[ex, 'sequence_3utr_start'] & pos<=these[ex, 'sequence_3utr_end']){
          utr[i] <- 3
        }
      }
      if(!is.na(these[ex, 'sequence_5utr_start'])){
        if(!is.na(these[ex, 'sequence_5utr_start'])){
          if(pos>=these[ex, 'sequence_5utr_start'] & pos<=these[ex, 'sequence_5utr_end']){
            utr[i] <- 5
          }
        }
      }
    } else{
                                        #intronic, where is it?
      start <- these[max(start.ok),]
      end <- these[min(end.ok),]
      if(start$rank>end$rank){
        intron[i] <- paste(end$'sequence_exon_stable_id',' (',end$rank,') - ',start$'sequence_exon_stable_id',' (',start$rank,')',sep="")
      }else{
        intron[i] <- paste(start$'sequence_exon_stable_id',' (',start$rank,') - ',end$'sequence_exon_stable_id',' (',end$rank,')',sep="")
      }
      
    }#else
    
  }#for

  return(data.frame(nearest, exon=exon, intron=intron, utr=utr))
}

#intra.ids[8] is exonic  



get.ensembl.gene.context <- function(nearest, ensmart){

  #change this when you've got some accessor methods
  nearest <- nearest@data
  exon  <- intron <- as.character(rep(NA,nrow(nearest)))

  #we'll set these to NA later.
  utr <- rep('',nrow(nearest))
  
  #only need context for intragenic features:
  intra.ids <- which( ( nearest$feat.GenomeMid > nearest$nearest.GenomeStart )
                      & ( nearest$feat.GenomeMid<nearest$nearest.GenomeEnd ) )

  
  #ok, grab those transcripts from ensembl.
  filters <- c('ensembl_gene_id')
  values <- nearest[intra.ids,"nearest.EnsemblGeneID"]

  attributes <- c(
                  'sequence_exon_stable_id',
                  'structure_exon_chrom_start',
                  'structure_exon_chrom_end',
                  'structure_exon_rank',
                  'sequence_exon_chrom_strand',
                  'sequence_3utr_start',
                  'sequence_3utr_end',
                  'sequence_5utr_start',
                  'sequence_5utr_end',
                  'ensembl_gene_id'
                 )

  
#  attributes<-c('sequence_exon_stable_id',
#                'sequence_exon_chrom_start',
#                'sequence_exon_chrom_end',
#                'rank',
#                'sequence_exon_chrom_strand',
#                '3_utr_start',
#                '3_utr_end',
#                '5_utr_start',
#                '5_utr_end',
#                'ensembl_gene_id',
#                )

  res <- getBM(attributes, filters=filters, values=values, mart=ensmart)

  n<-0
  len<-length(intra.ids)
  

  for(i in intra.ids){
    cat(i,"\n")
    n<-n+1
    cat(n,'of',len,"\n" )
    ts<-nearest$nearest.EnsemblGeneID[i]
    these <- res[res$'ensembl_gene_id'== ts,]
    pos <- nearest$feat.GenomeMid[i]
    start.ok <- which(these$'exon_chrom_start'<=pos)
    end.ok <- which(these$'exon_chrom_end'>= pos)
    ex <- intersect(start.ok, end.ok)
    strand <- these$'exon_chrom_strand'[1]
    
    if(length(ex)!=0){
      
      #exonic, list all exons in which it appears
      exon[i] <- paste(these[ex,'exon_stable_id'],' (' ,these[ex,'rank'], ')',sep="", collapse=", ")
                                        #see if it's in the utr
      these.ex <- these[ex,]
      for (r in 1:nrow(these.ex)){
        this.utr <- 'NA'
        if(!is.na(these.ex[r,'3utr_start']) &
           pos>=these.ex[r, '3utr_start'] &
           pos<=these.ex[r, '3utr_end']){
          this.utr <- '3'
        }
        if(!is.na(these.ex[r, '5utr_start']) &
           pos>=these.ex[r, '5utr_start'] &
           pos<=these.ex[r, '5utr_end']){
          this.utr <- '5'
        }
        
        if(utr[i]==""){utr[i] <- this.utr} else {utr[i] <- paste(utr[i],this.utr,sep=", ")}
        
    }#for

      
  } else{

    #this is going to be more complicated than I thought
    #because what we get back is the total set of exons from
    #all transcripts of the gene (with dups removed)
    #the ranks of the exons may differ between transcripts
    #just make a note of whether it's in the first intron
    #for now and fix properly later.
    
    start <- these[max(start.ok),]
    end <- these[min(end.ok),]
    
    intron[i]<-"INTRONIC"
    if(end$rank %in% c(1,2) & start$rank %in% c(1,2)){
     intron[i] <- 'FIRST INTRON'
   }
    
  }#else
    
  }#for

  #fix the utr list for non-exonic ones
  utr[which(utr=='')] <- NA
  
  return(data.frame(nearest, exon=exon, intron=intron, utr=utr))
}
