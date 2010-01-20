###
# R script to generate a graphical representation of
# the interactions between researchers at the IoP
# See Rakefile.iop_connections for more details

# CassJ <cassjohnston@gmail.com> Jan 2010.

library(igraph)

# load the connectivity data
data <- read.csv("Sheet1.csv", header=T, stringsAsFactors=F)
rownames(data) <- data[,1]
data <- data[,-1]


# change 'x' to 1 and NA or '' to 0
# there's probably a more Rish way of doing this.
for (i in 1:nrow(data)){
  inds <- union(which(is.na(data[i,])), which(data[i,]==""))
  data[i,]<-1
  data[i, inds]<-0
}


#sort the metadata into affiliations
metadata <- metadata[order(metadata[,"Affiliation"]),]

#sort the data in the same order
inds <- rownames(metadata)
data<-data[inds,]
data<-data[,inds]

# load the metadata about people 
metadata <- read.csv("Sheet2.csv", header=T, stringsAsFactors=F, row.names=1)
people <- rownames(data)
affiliations <- unique(metadata[,"Affiliation"])
cols <- rainbow(length(affiliations))
names(cols) <- affiliations

#create an empty graph
iop.graph <- graph.formula()


#populate it with vertices
for (i in people){
 col <- cols[metadata[i,"Affiliation"]]
 iop.graph<-add.vertices(iop.graph, 1, label=i, color=col)
}

#Note that Vertex and Edge ids are 0-based

#And draw edges between people
for(i in 1:length(people)){
  for(j in 1:length(people)){
       if(data[i,j]==1){
       iop.graph <- add.edges(iop.graph, c((i-1),(j-1)))
       }
 }
}


iop.graph <- set.graph.attribute(iop.graph, "layout", layout.circle(iop.graph))


library(Cairo)
Cairo(file="collaboration.pdf", type="pdf", width=11.7, height=8.3, title="Collaboration in the IoP", units="in", dpi=300)
plot(iop.graph)

title(main="Collaboration within the IoP", cex.main=2)

savefont <- par(font=3) 
legend("bottomleft", legend=names(cols), fill=cols, bty="n", cex=0.7)


namemap <- paste(rownames(metadata), "     ", metadata[,1], metadata[,2])
legend("bottomright", legend=namemap, bty="n", cex=0.7)
par(savefont)

dev.off()
