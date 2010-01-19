###
# R script to generate a graphical representation of
# the interactions between researchers at the IoP
# See Rakefile.iop_connections for more details


data <- read.csv("Sheet1.csv", header=T)
rownames(data) <- data[,1]
data <- data[,-1]



