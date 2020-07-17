args <- commandArgs(TRUE)
rf <- as.character(args[1])
ABC <- as.character(args[2])
outname <- as.character(args[3])

library(GeneArchEst)

rf <- readRDS(rf)
ABC <- readRDS(ABC)

scale <- hyperparameter_regression_on_ABC(ABC$ABC_res, rf$pi$parameter_density$pi, acceptance_threshold = 0.0005)

saveRDS(scale, outname)