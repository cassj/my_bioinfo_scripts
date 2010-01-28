source("scripts/qw.R")
source("scripts/liftOver.R")

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))

data <- read.csv(filename, header=F, sep="\t")

to.map<-data[,1:3]
colnames(to.map) <- c("chr","start","end")
mapped <- liftOver(to.map , chain.file="lib/mm8ToMm9.over.chain")
data[,1:3] <- mapped[,qw(chr,start,end)] 

outfile=sub('mm8','mm9',filename)

write.table(data, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
