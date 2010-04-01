#call like R --vanilla --args filename=\"thing\" start=1 end=10 outfile=\"otherthing\"
args<-commandArgs()
eval(parse(text=args[grep('filename', args)]))
eval(parse(text=args[grep('outfile', args)]))

eval(parse(text=args[grep('start', args)]))
eval(parse(text=args[grep('end', args)]))

data <- read.csv(filename, as.is=T)

#delete any existing version
unlink(outfile)

for (i in start:end){
  write(paste('>', data$names[i], "\n", data$sequence[i], sep=""),append=T, file=outfile)
}


