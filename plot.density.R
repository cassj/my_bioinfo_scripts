
plot.mat.density<-function(x,log=F){
 if(log){x<-log2(x)}
 plot(density(x[,1]), main="Density Distributions")
 for(i in 2:ncol(x)){lines(density(x[,i]),col=i)}
 legend("topright", legend=colnames(x), fill=1:ncol(x))
}

