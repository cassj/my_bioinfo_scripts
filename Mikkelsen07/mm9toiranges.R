source("scripts/qw.R")
library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))

data <- read.csv(filename, header=F, sep="\t")

colnames(data) <- qw(Chr, Start, End, Name, Score)

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
data[,"Name"]<-paste(paste(data[,"Chr"],data[,"Start"],sep=":"),data[,"End"], sep="-")

#some inexpliably have a row of NA at the bottom.
data<-data[!is.na(data[,"Chr"]),]

rd <- RangedData(ranges = IRanges(
                   start= data$Start,
                   end = data$End,
                   ),
                 names = as.character(data$Name),
                 space = as.character(data$Chr),
                 values = data[,"Score"]
                 )

outfile=sub('.mm9.bed','.RangedData.R',filename)

#And save the result
save(rd, file=outfile)


