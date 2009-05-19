

setClassUnion("matrixy", c("matrix", "data.frame"))

#Virtual base class
setClass("featureComp",
         representation(
                        genomeBuild="character",
                        species="character",
                        genomeFeature="character",
                        compared.to="character",
                        feature.context="numeric",
                        just.data.2="logical",
                        notes="character",
                        "VIRTUAL"
                        )

         )


setClass("nearestFeature",
         representation(
                        data="matrixy",
                        field.1 = "character",
                        field.2 = "character",
                        n="numeric"
                        ),
         contains="featureComp"
         )


setClass("overlappingFeature",
         representation(
                        data="matrixy",
                        up="numeric",
                        down="numeric"
                        ),
         contains="featureComp"
         )



# user-friendly constructors

#basic constructor, just passes all args to new
nF <- nearestFeature  <- function(...){
  new ("nearestFeature", ...)
}
oF <- overlappingFeature  <- function(...){
  new ("overlappingFeature", ...)
}
