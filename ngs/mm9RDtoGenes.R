library(IRanges)

#call like R --vanilla --args filename=\"thing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
load(filename)

library(ChIPpeakAnno)

#load mouse transcripts
data(TSS.mouse.NCBIM37)

#need to alter the names so they match ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

data.annot <- annotatePeakInBatch(rd, AnnotationData=TSS.mouse.NCBIM37)

newfile<-sub("RangedData", "AnnoRangedData", filename)
save(data.annot, file=newfile)
