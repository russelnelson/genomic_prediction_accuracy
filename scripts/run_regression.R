args <- commandArgs(TRUE)
rf <- as.character(args[1])
ABCfile <- as.character(args[2])
outname <- as.character(args[3])

library(GeneArchEst)

rf <- readRDS(rf)
nsites <- rf$sites$parameter_density$sites
rm(rf)
gc();gc()

res <- data.table::fread(ABCfile, header = F)
colnames(res) <- c("sites", "d.f", "scale", "h", GeneArchEst::names_diff_stats)

scale <- hyperparameter_regression_on_ABC(res, nsites, acceptance_threshold = 0.005, formula = scale ~ sites,
                                          parameter_transforms = reasonable_transform(c("scale", "sites"))$forward,
                                          parameter_back_transforms = reasonable_transform(c("scale", "sites"))$back)

saveRDS(scale, outname)
