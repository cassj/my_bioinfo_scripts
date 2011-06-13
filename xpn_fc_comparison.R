#!/usr/local/bin/Rscript


neurons.file <- "../mla_neurons_lumixpn/results/limma_results.csv"
astrocytes.file <- '../ns5_vs_ns5dastro_lumixpn/results/limma_rd.csv'

neurons<- read.csv(neurons.file)
astrocytes <- read.csv(astrocytes.file)

ensid.a<-astrocytes[,"ensembl.gene.id"]
ensid.n<-neurons[,"EnsemblID"]

common.ens<-intersect(ensid.n, ensid.a)

astro.common<-astrocytes[ensid.a %in% common.ens,]
neurons.common<-neurons[ensid.n %in% common.ens,]

#order by absolute value of logFC

astro.common<-astro.common[order(abs(astro.common[,6]), decreasing=TRUE),]
neurons.common<-neurons.common[order(abs(neurons.common[,10]), decreasing=TRUE),]

#remove duplicates
astronodup<-astro.common[!duplicated(astro.common[,28]),]
neuronsnodup<-neurons.common[!duplicated(neurons.common[,2]),]



#change row names to ensembl id so that it is ordered by that on the pretty plot

rownames(astronodup) <- astronodup[,28]
rownames(neuronsnodup) <- neuronsnodup[,2]
all.ids<-rownames(astronodup)

library(biomaRt)
ensmart<- useMart("ensembl",dataset="mmusculus_gene_ensembl")
results<- getBM(filters="ensembl_gene_id", values=rownames(astronodup), attributes=c("ensembl_gene_id", "go_biological_process_id"), mart=ensmart)

results.go<- unique(results[results[,"go_biological_process_id"]=="GO:0022008","ensembl_gene_id"])
#transcription factor 0003700
#cell cycle 0007049
#neurogenesis 0022008
#regulation of neurogenesis 0050767
#positive regulation of neurogenesis 0050769
#negative regulation of neurogenesis 0050768
#neuronal migration 0001764

postscript(file="results/dotplotNSvsAstrotranscriptionfactors.eps", paper="special", height=10, width=10)

# for plotting just TFs
all.ids<-intersect(all.ids,results.go)

plot(astronodup[all.ids,"logFC"], neuronsnodup[all.ids,"logFC"], xlab="Fold Change between NS and Astrocytes", ylab="Fold Change between NS and Neurons", pch=".")


# subset data for plotting significant points
sigastro.ids <- rownames(astronodup[astronodup[,10] <= 0.0001, ])
signeurons.ids <- rownames(neuronsnodup[neuronsnodup[,14] <=0.0001, ])
sig.both.ids <- intersect(sigastro.ids, signeurons.ids)

#For plotting just TFs
sigastro.ids<-intersect(sigastro.ids, results.go)
signeurons.ids<-intersect(signeurons.ids, results.go)
sig.both.ids<-intersect(sig.both.ids, results.go)

# colour significant astrocyte points in red
points(astronodup[sigastro.ids,"logFC"], neuronsnodup[sigastro.ids,"logFC"],  col="red")

# colour significant neuron points in blue
points(astronodup[signeurons.ids,"logFC"], neuronsnodup[signeurons.ids,"logFC"], col="blue")

# colour sig in both points purple
points(astronodup[sig.both.ids,"logFC"], neuronsnodup[sig.both.ids,"logFC"], col="magenta")

#plot go transcription factors in green
#points(astronodup[results.go,"logFC"], neuronsnodup[results.go,"logFC"], pch="x", col="green")

# split the graph into quarters
xrange <- range(astronodup[,"logFC"])
yrange <- range(neuronsnodup[,"logFC"])
lines(xrange,c(0,0), lty=2)
lines(c(0,0), yrange, lty=2)

#add Ensembl IDs to points when you click on the plot. Obviously this won't work with poscript output
#will investigate a better way of doing this
identify(astronodup[all.ids, "logFC"], neuronsnodup[all.ids,"logFC"], labels=neuronsnodup[all.ids,"symbol"])


dev.off()

# subset data for plotting significant points but less harsh
sigastro.ids <- rownames(astronodup[astronodup[,10] <= 0.01, ])
signeurons.ids <- rownames(neuronsnodup[neuronsnodup[,14] <=0.01, ])
sig.both.ids <- intersect(sigastro.ids, signeurons.ids)


astro.up.ids <- rownames(astronodup[astronodup[,"logFC"]>=1,])
astro.down.ids <- rownames(astronodup[astronodup[,"logFC"]<=-1,])

neuron.up.ids <- rownames(neuronsnodup[neuronsnodup[,"logFC"]>=1,])
neuron.down.ids <- rownames(neuronsnodup[neuronsnodup[,"logFC"]<=-1,])

#get genes that are not changing
astro.nochange.ids <-rownames(astronodup[(astronodup[,"logFC"]>=-1 & astronodup[,"logFC"]<=1),])
#neurons.nochange.ids <-rownames(neuronsnodup[(neuronsnodup[,"logFC"]>=-1 & neuronsnodup[,"logFC"]<=1),])
neurons.nochange.ids <- rownames(neuronsnodup[abs(neuronsnodup[,"logFC"])<=1,])

astro.sig.up<-astronodup[intersect(astro.up.ids, sigastro.ids),]
astro.sig.up<-astro.sig.up[order(astro.sig.up[,"logFC"], decreasing=T),]
write.csv(astro.sig.up, file="results/astro_sig_up.csv")

astro.sig.down<-astronodup[intersect(astro.down.ids, sigastro.ids),]
astro.sig.down<-astro.sig.up[order(astro.sig.down[,"logFC"]),]
write.csv(astro.sig.down, file="results/astro_sig_down.csv")

neuron.sig.up<-neuronsnodup[intersect(neuron.up.ids, signeurons.ids),]
neuron.sig.up<-neuron.sig.up[order(neuron.sig.up[,"logFC"], decreasing=T),]
write.csv(neuron.sig.up, file="results/neuron_sig_up.csv")

neuron.sig.down<-neuronsnodup[intersect(neuron.down.ids, signeurons.ids),]
neuron.sig.down<-neuron.sig.up[order(neuron.sig.down[,"logFC"]),]
write.csv(neuron.sig.down, file="results/neuron_sig_down.csv")

#get genes that are only changing in astro or in neurons
neuron.unique.changing.up<-neuronsnodup[intersect(intersect(neuron.up.ids, astro.nochange.ids), signeurons.ids),]
neuron.unique.changing.down<-neuronsnodup[intersect(intersect(neuron.down.ids, astro.nochange.ids), signeurons.ids),]

astro.unique.changing.up<-astronodup[intersect(intersect(astro.up.ids, neurons.nochange.ids), sigastro.ids),]
astro.unique.changing.down<-astronodup[intersect(intersect(astro.down.ids, neurons.nochange.ids), sigastro.ids),]

#order them by logFC
neuron.unique.changing.up<-neuron.unique.changing.up[order(neuron.unique.changing.up[,"logFC"], decreasing=T),]
neuron.unique.changing.down<-neuron.unique.changing.down[order(neuron.unique.changing.down[,"logFC"], decreasing=T),]

astro.unique.changing.up<-astro.unique.changing.up[order(astro.unique.changing.up[,"logFC"], decreasing=T),]
astro.unique.changing.down<-astro.unique.changing.down[order(astro.unique.changing.down[,"logFC"], decreasing=T),]

#print them out
write.csv(neuron.unique.changing.up, file="results/neuron.unique.changing.up.csv")
write.csv(neuron.unique.changing.down, file="results/neuron.unique.changing.down.csv")

write.csv(astro.unique.changing.up, file="results/astro.unique.changing.up.csv")
write.csv(astro.unique.changing.down, file="results/astro.unique.changing.down.csv")

#find genes that change in both - either up or down
both.up.ids  <- intersect(sig.both.ids, intersect(neuron.up.ids, astro.up.ids))
both.down.ids <- intersect(sig.both.ids, intersect(neuron.down.ids, astro.down.ids))
up.astro.down.neuron.ids <- intersect(sig.both.ids, intersect(astro.up.ids, neuron.down.ids))
down.astro.up.neuron.ids <- intersect(sig.both.ids, intersect(astro.down.ids, neuron.up.ids))
both.opp.ids<-c(up.astro.down.neuron.ids, down.astro.up.neuron.ids)

colnames(astronodup)<-paste("astro", colnames(astronodup),sep="_")
colnames(neuronsnodup)<-paste("neuron", colnames(neuronsnodup), sep="_")
both<-cbind(astronodup[all.ids,], neuronsnodup[all.ids,])

both<-both[order(both[,"neuron_logFC"], both[,"astro_logFC"], decreasing=T ),]
write.csv(both[both.up.ids,], file="results/both_sig_up.csv")

both<-both[order(both[,"neuron_logFC"], both[,"astro_logFC"], decreasing=F ),]
write.csv(both[both.down.ids,], file="results/both_sig_down.csv")

write.csv(both[both.opp.ids,], file="results/both_sig_opp.csv")


#save absolute values

#astro.changing<-astronodup.p[abs(astronodup.p[,6]) >= 1, ]
#neurons.changing<-neuronsnodup.p[abs(neuronsnodup.p[,10])>=1, ]



#common<-intersect(astro.changing[,28], neurons.changing[,2])


#how many genes are up and down regulated in neurons and/or astrocytes
#astro.changing.up<-astro.changing[astro.changing[,6] >= 1, ]
#neurons.changing.up<-neurons.changing[neurons.changing[,10]>=1, ]

#astro.changing.down<-astro.changing[astro.changing[,6] <= 1, ]
#neurons.changing.down<-neurons.changing[neurons.changing[,10] <=1, ]

#common.up<-intersect(astro.changing.up[,28], neurons.changing.up[,2])
#common.down<-intersect(astro.changing.down[,28], neurons.changing.down[,2])


#make a csv file of genes that are changing in both neurons and astrocytes
#common.genes.neurons <- neurons.changing[common,]
#common.genes.astro <- astro.changing[common,]

#common.genes.dat <- cbind(rownames=common,common.genes.neurons[,10], common.genes.neurons[,14], common.genes.astro[,6], common.genes.astro[,10])

#write.csv(common.genes.dat, file="commongeneschange.csv")


#save commonly significant genes
#astro.common<-astronodup[common,]
#neurons.common<-neuronsnodup[common,]

#trying to get a list of genes that change in both astrocytes and neurons with a fold change column.
#common<-intersect(astro.common[,28], neurons.common[,2])

#save changing genes
#write.csv(astro.changing, file="results/astro.changing.csv")
#write.csv(neurons.changing, file="results/neurons.changing.csv")

#save common genes
#write.csv(astro.common, file="results/astro.common.csv")
#write.csv(neurons.common, file="results/neurons.common.csv")







