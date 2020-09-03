library(GeneArchEst)

args <- commandArgs(TRUE)
x <- as.character(args[1])
out <- as.character(args[2])

x <- data.table::fread(x, header = F)
x <- as.data.frame(x)
colnames(x) <- c("chr", "position")
windows <- mark_windows(x, 50, "chr")
saveRDS(windows, out)
