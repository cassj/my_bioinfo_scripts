#build a new PWM from our data


library(rGADEM)
library(BSgenome.Mmusculus.UCSC.mm9)

eval(parse(text=args[grep('filename', args)]))


#filename<-"Macs/run4/top10K.fa"
