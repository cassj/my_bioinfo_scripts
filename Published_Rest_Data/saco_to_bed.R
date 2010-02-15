#call like R --vanilla --args filename=\"thing\"

source("scripts/qw.R")
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
data<- read.table(filename,sep="\t", header=T)

bed <- data[,qw(Chr,Start,End)]
newfile <- 'results/saco_mm9.bed'
write.table(bed, file=newfile, sep="\t", col.names=F, row.names=F, quote=F)
