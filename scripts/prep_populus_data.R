# read and prep the populus data
## read
library(GeneArchEst)
args <- commandArgs(TRUE)
dat <- as.character(args[1])
phenos <- as.character(args[2])

# read in phenotypes, figure out which have data for our phenotype
phenos <- data.table::fread(phenos)



## clean and save
