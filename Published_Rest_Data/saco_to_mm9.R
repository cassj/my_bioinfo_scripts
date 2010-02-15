source('scripts/qw.R')
source('scripts/liftOver.R')

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data<- read.table(filename,sep="\t", header=T)

#All of the tags are 16 bases long.
#I have to assume the positions given are 1-based
to.map<-data.frame(chr=as.character(data$Chr),
                   start=as.numeric(data$Start-1))
to.map <- cbind(to.map, end=to.map$start+16)

#mm5 to mm8, then mm8 to mm9 cos there isn't a direct chain file.
mapped.mm8 <- liftOver(to.map, chain.file="lib/mm5ToMm8.over.chain")
mapped.mm9 <- liftOver(mapped.mm8, chain.file="lib/mm8ToMm9.over.chain")

#and stick the results back in the dataframe
data[,"Start"] <- mapped.mm9[,"start"]
data[,"Chr"] <- mapped.mm9[,"chr"]
data <- data.frame(data, End=mapped.mm9[,"end"])

#liftover is only approximate, so the actual sequence match may be
#a few bases away from the liftover position.


#save the results
outfile=sub('html','txt',filename)
write.table(data, file=outfile, sep="\t", col.names=F, row.names=T, quote=F)
