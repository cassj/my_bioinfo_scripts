#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))

data <- read.csv(filename, header=F, sep="\t")

#the only reason this isn't bed is that it is missing a name col
data<-cbind(data[,1:3], NA, data[,4])

outfile=sub('.txt','.mm8.bed', sub('publication','results', filename))

write.table(data, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
