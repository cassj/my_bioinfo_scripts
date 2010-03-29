library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data.annot <- get(load(filename))

data<-as.data.frame(data.annot)
colnames(data)[1]<-"chr"
colnames(data) <- gsub('values.', '', colnames(data))
data <- data[order(data$neg10log10pVal, decreasing=T),]
pval <- 10^(data$neg10log10pVal/-10)
data <- cbind(data, pval)
  
newfile<-sub("\\.R", ".csv", filename)
write.csv(data, file=newfile, row.names=F)
