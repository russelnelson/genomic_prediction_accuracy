args <- commandArgs(TRUE)
sim_dat <- as.character(args[1])
geno_pheno <- as.character(args[2])
outname <- as.character(args[3])
data_type <- as.character(args[4])

library(GeneArchEst)



sim_dat <- readr::read_delim(sim_dat,
                             "\t", escape_double = FALSE, col_names = FALSE,
                             trim_ws = TRUE)
colnames(sim_dat) <- c("pi", "scale", "d.f", "h", names_descriptive_stats)


x <- readRDS(geno_pheno)
phenos <- x$phenos$p
meta <- x$meta
if(data_type == "imputed"){
  x <- x[[3]]
}
if(data_type == "full"){
  x <- x[[1]]
}

# note, not naming and removing the iter column, and not checking for that error, produces an uninformative error
# in the quantforresterror. Might gather some example data and post it as an issue.
pred <- hyperparameter_random_forest(x = x, meta = meta, phenos = phenos, num_trees = 100000,
                                    sims = sim_dat, num_threads = 24, hyperparameter_to_estimate = "pi")


saveRDS(pred, outname)
