args <- commandArgs(TRUE)
sim_dat <- as.character(args[1])
metafile <- as.character(args[2])
windowfile <- as.character(args[3])
gmmatfile <- as.character(args[4])
grmfile <- as.character(args[5])
phenofile <- as.character(args[6])
outname <- as.character(args[7])

library(GeneArchEst)


# read in the sim GWAS summaries
sim_dat <- readr::read_delim(sim_dat,
                             "\t", escape_double = FALSE, col_names = FALSE,
                             trim_ws = TRUE)
colnames(sim_dat) <- c("pi", "scale", "d.f", "h", names_descriptive_stats)


# read in the windows
pass_windows <- readRDS(windowfile)

# read in the G matrix
pass_G <- data.table::fread(grmfile, header = F)
pass_G <- as.matrix(pass_G)

# read in the subset metadata
meta <- as.data.frame(data.table::fread(metafile, header = F))
colnames(meta) <- c("chr", "position")

# read in the phenotypes
phenos <- readRDS(phenofile)$phenos$p


# note, not naming and removing the iter column, and not checking for that error, produces an uninformative error
# in the quantforresterror. Might gather some example data and post it as an issue.
pred <- hyperparameter_random_forest(x = NULL, meta = meta, phenos = phenos, num_trees = 100000,
                                    sims = sim_dat, num_threads = 24, hyperparameter_to_estimate = "pi", 
                                    pass_windows = pass_windows, pass_G = pass_G, GMMAT_infile = gmmatfile)


saveRDS(pred, outname)
