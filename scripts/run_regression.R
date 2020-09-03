args <- commandArgs(TRUE)
rf <- as.character(args[1])
ABCfile <- as.character(args[2])
outname <- as.character(args[3])

library(GeneArchEst)

rf <- readRDS(rf)
pi <- rf$pi$parameter_density$pi
rm(rf)
gc();gc()

res <- data.table::fread(ABCfile, header = F)
colnames(res) <- c("pi", "d.f", "scale", "h", GeneArchEst::names_diff_stats)

scale <- hyperparameter_regression_on_ABC(res, pi, acceptance_threshold = 0.0005)

saveRDS(scale, outname)
