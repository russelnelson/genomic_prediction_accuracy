# read and prep the populus data
args <- commandArgs(TRUE)
dat <- as.character(args[1])
#dat <- "../genomic_prediction_accuracy/data/populus_data/cp_test.vcf"
phenos <- as.character(args[2])
#phenos <- "../genomic_prediction_accuracy/data/populus_data/sorted_phenotypes.txt"

# read in phenotypes, figure out which have data for our phenotype
phenos <- data.table::fread(phenos)
goods <- which(!is.na(phenos$`Stomatal density (SD)`))
good.samps <- phenos$JGIname[goods]
phenos <- phenos$`Stomatal density (SD)`[goods]

#==================format, clean, and save genotypes as a FBM====================
## find the cols to import
systype <- Sys.info()[1]
if(systype == "Windows"){
  filemeta <- shell(paste("grep '##'", dat,  "| wc -l"), intern = T)
}
if(!exists("filemeta")) {
  filemeta <- system(paste("grep '##'", dat,  "| wc -l"), intern = T)
}
filemeta <- as.numeric(filemeta)
header <- unlist(data.table::fread(dat, skip = filemeta, nrows = 1, header = F))
keep.cols <- which(header %in% good.samps)
## read in
genos <- bigstatsr::big_read(dat, select = keep.cols, skip = filemeta + 1, type = "integer")
genos <- genos$save()

#==================get and save metadata===============
meta <- data.table::fread(dat, select = 1:2, skip = filemeta + 1)
colnames(meta) <- c("chr", "position")
saveRDS(list(meta = meta, phenos = list(p = phenos)), paste0(dat, ".meta.RDS"))
