source("scripts/qw.R")

# Kee Yew has already filtered the data on
# detection pvalue and background corrected it.
# Bugger.

options(scipen=10)

#load data
data <- read.csv("expression_data/Astro-REST-DN-raw.csv", as.is=TRUE)

#grab expression value cols
cols <- qw(dn, dn, dn, dn, ev, ev, ev, ev)
rownames(data)<-data[,1]
n<-length(cols)
data<-data[,2:(n+1)]
ev<-which(cols=="ev")
dn<-which(cols=="dn")


#set values <10 to 10
fix10<-function(x){
  x[x<10]<-10
  return(x)
}

data<-apply(data,2,fix10)

#quantile normalise. The normaliseIllumina method seems to cope with 
#an expressionset, which save the hassle of trying to build a valid
#BSData object
library(affy)
library(beadarray)

es <- new("ExpressionSet", exprs = data)
E = normaliseIllumina(es, method="quantile", transform="log2")
data <- exprs(E)


library(limma)

design<-matrix(0,nrow=(ncol(data)), ncol=2)
colnames(design)<-c("ev","dn")
design[ev,"ev"]<-1
design[dn,"dn"]<-1


fit<-lmFit(data, design)

cont.matrix<-makeContrasts(dnvsev=dn-ev, levels=design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)
write.fit(ebFit, file="expression_data/limma_results.csv", adjust="BH")
data<-read.table("expression_data/limma_results.csv", sep="\t", header=T)

data<- topTable(ebFit, number=nrow(data))
write.csv(data,"expression_data/limma_results.csv")




