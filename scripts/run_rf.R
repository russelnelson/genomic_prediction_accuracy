args <- commandArgs(TRUE)
sim_dat <- as.character(args[1])
geno_pheno <- as.character(args[2])
outname <- as.character(args[3])
peak_delta <- as.numeric(args[4])
peak_pcut <- as.numeric(args[5])

library(GeneArchEst)



sim_dat <- readr::read_delim(sim_dat,
                             "\t", escape_double = FALSE, col_names = FALSE,
                             trim_ws = TRUE)
colnames(sim_dat) <- c("sites", "scale", "d.f", "h", names_descriptive_stats)


x <- readRDS(geno_pheno)
phenos <- x$phenos
meta <- x$meta
meta <- as.data.frame(meta)
x <- x$genos

pred <- hyperparameter_random_forest(x = x, meta = meta, phenos = phenos$p, num_trees = 10000, importance = "none", hyperparameter_to_estimate = "sites",
                                    peak_delta = peak_delta, peak_pcut = peak_pcut, center = TRUE,
                                    sims = sim_dat, num_threads = 24, phased = FALSE)


saveRDS(pred, outname)
