#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data<- read.csv(filename)

pVal<-10^(data[,"values.neg10log10pVal"]/-10)
ord<-order(pVal)
data<-cbind(data, pVal)
data<-data[ord,]



#save the results
outfile=sub('RangedData','RangedDatapVal',filename)
write.csv(data, file=outfile, row.names=F)
