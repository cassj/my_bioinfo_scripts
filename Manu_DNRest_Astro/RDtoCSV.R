library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data.annot <- get(load(filename))

data<-as.data.frame(data.annot)
colnames(data)[1]<-"chr"

newfile<-sub("\\.R", ".csv", filename)
write.csv(data, file=newfile, row.names=F)
