library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))

library(ChIPpeakAnno)

#load mouse transcripts
data(TSS.mouse.NCBIM37)

#
