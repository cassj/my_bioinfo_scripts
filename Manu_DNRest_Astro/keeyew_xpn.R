source("scripts/qw.R")

# Following Kee Yew's protocol
# Data has already been background corrected.
# I don't know how, check with KY, but it has -ve values

options(scipen=10)

#load data
data <- read.csv("expression_data/Astro-REST-DN-raw.csv", as.is=TRUE)

cols <- qw(dn, dn, dn, dn, ev, ev, ev, ev)
rownames(data)<-data[,1]
n<-length(cols)
data<-data[,2:(n+1)]

#set values <10 to 10
fix10<-function(x){
  x[x<10]<-10
  return(x)
}

data<-apply(data,2,fix10)

#According to the GEO GPL4234 annotation for the Sentrix ref6 array
#there are 47,000 transcripts represented on the array.
#for the esc data we have 20,433. hmm. shite.
#for nsc 21297
#for astro 15060
#so KeeYew has already filtered on detection pvalue 



# normalise per chip to median 

chip.meds<-apply(data, 2, median)
data<-data/rep(chip.meds, each=nrow(data))

# normalise per gene to median

gene.meds<-apply(data, 1, median)
data<-data/gene.meds
save(data, file="expression_data/normalised_data.RData")


###
#de: just a 2-sided t-test with Welch's estimate 

library(stats)
dn<-which(cols=="dn")
ev<-which(cols=="ev")


fc<-apply(data, 1, function(x){mean(x[dn])/mean(x[ev])})
names(fc)<-rownames(data)

tests<-apply(data,1,function(x){ t.test(x[dn],x[ev] )  })
save(tests, file="expression_data/ky_tests.RData")

t.stats<-unlist(lapply(tests,function(x){x$statistic}))
names(t.stats)<-rownames(data)
t.stats<-sort(t.stats, decreasing=TRUE)
save(t.stats, file="expression_data/ky.tstats.RData")

p.vals<-unlist(lapply(tests,function(x){x$p.value}))
names(p.vals)<-rownames(data)
p.vals<-sort(p.vals)
save(p.vals, file="expression_data/ky_pvals.RData")
ky.ord<-names(p.vals)




###
# Benjamini-Hochberg FDRs

library(multtest)

bh.fdr<-mt.rawp2adjp(p.vals, "BH")[[1]]
rownames(bh.fdr)<-names(p.vals)
bh.fdr<-cbind(FC=fc[rownames(bh.fdr)], bh.fdr)
bh.fdr<-cbind(LogFC=log2(fc[rownames(bh.fdr)]), bh.fdr)

save(bh.fdr, file="expression_data/ky_results.RData")

kyfdr.05<-rownames(bh.fdr[which(bh.fdr[,"BH"]<0.05),])
write.csv(bh.fdr, file="expression_data/ky_results.csv")


