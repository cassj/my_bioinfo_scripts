#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))

#the only reason this isn't bed is that it has colnames
#read it in, don't write it out
data <- read.csv(filename, header=T, sep="\t")

#and for consistency, we want the chr names to be in the form chrN
data[,1]<-paste('chr',data[,1], sep="")

#and the same num of cols as the windowinterval data
data<-cbind(data,name=NULL, score=NULL)

outfile=sub('.txt','.mm8.bed', sub('publication','results', filename))
write.table(data, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
