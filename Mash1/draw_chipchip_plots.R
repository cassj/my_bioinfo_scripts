#quick script to generate some plots for DC

#"Norm IP Signal_0" vs Norm WCE Signal_0. (cols N and O)
# Red for probes in bound regions, information that is contained in column "V" - "Is In BoundRegion"


telenc.dir <- "/home/cassj/work/Telencephalon final"
ns.dir <- "/home/cassj/work/NS5 final"

telenc.files <- paste(telenc.dir, grep("ChIP_Analytics_1.3.1_Output_Probe",dir(telenc.dir), value=T), sep="/")
ns.files <- paste(ns.dir, grep("ChIP_Analytics_1.3.1_Output_Probe",dir(ns.dir), value=T), sep="/") 


do.plot<-function(f){
  for (i in 1:length(f)){
    cat("plotting file ",i, "\n")
    this.data<-read.table(f[i], sep="\t", header=T)
    x <- this.data[,"Norm.WCE.Signal_0"]
    y <- this.data[,"Norm.IP.Signal_0"]
    
    #a few have NaN values. Don't know why, but we'll take them out.
    rm.inds <- -1 * union(which(is.na(x)), which(is.na(y)))
    if(length(rm.inds)>0){
      x <- x[rm.inds]
      y <- y[rm.inds]
      this.data <- this.data[rm.inds,]
    }
    
    bound <- which(this.data[,"Is.In.BoundRegion"]==1)
    not.bound <-  which(this.data[,"Is.In.BoundRegion"]==0)
    
    x<-log2(x)
    y<-log2(y)

    #plot the axes to begin with
    if(i==1){
      plot(c(0,16),c(0,16), col="white", xlab="", ylab="")
    }
    
    points(x[not.bound], y[not.bound], pch=".", )
    points(x[bound],y[bound],pch=".", col="red")
    
  }
       
}


postscript(file=paste(telenc.dir, "graph.ps", sep="/"))
do.plot(telenc.files)
dev.off()

postscript(file=paste(ns.dir, "graph.ps", sep="/"))
do.plot(ns.files)
dev.off()
