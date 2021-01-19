# read and prep the populus data
args <- commandArgs(TRUE)
dat <- as.character(args[1])
phenos <- as.character(args[2])
phenos <- read.table(phenos)
phenos <- phenos[,1]

#==================format, clean, and save genotypes as a FBM====================
# find the cols to import
systype <- Sys.info()[1]
if(systype == "Windows"){
  filemeta <- shell(paste("grep '##'", dat,  "| wc -l"), intern = T)
}
if(!exists("filemeta")) {
  filemeta <- system(paste("grep '##'", dat,  "| wc -l"), intern = T)
}
filemeta <- as.numeric(filemeta)
header <- unlist(data.table::fread(dat, skip = filemeta, nrows = 1, header = F))

cat("reading in input:\n", dat, "\n")

# read in
genos <- bigstatsr::big_read(dat, select = 10:length(header), skip = filemeta + 1, type = "integer")
genos <- genos$save()

#==================get and save metadata===============
cat("reading in metadata\n")
meta <- data.table::fread(dat, select = 1:2, skip = filemeta + 1)
colnames(meta) <- c("chr", "position")
saveRDS(list(meta = meta, phenos = list(p = phenos)), paste0(dat, ".meta.RDS"))
writeLines(header[-c(1:9)], paste0(dat, ".keep.samples.txt"), sep = "\n")
