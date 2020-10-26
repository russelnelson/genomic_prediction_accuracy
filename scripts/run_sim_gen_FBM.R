library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
df_func <- function(x) runif(x, 1, 100)


args <- commandArgs(TRUE)
genofile <- as.character(args[1])
ABCfile <- as.character(args[2])
metafile <- as.character(args[3])
windowfile <- as.character(args[4])
gmmatfile <- as.character(args[5])
grmfile <- as.character(args[6])
hmean <- as.numeric(args[7])
hsd <- as.numeric(args[8])
outname <- as.character(args[9])

iters <- 1

# read in the genotypes
x <- bigstatsr::big_attach(genofile)


# read in the ABC results
res <- data.table::fread(ABCfile, header = F)
colnames(res) <- c("sites", "d.f", "scale", "h", GeneArchEst::names_diff_stats)
res <- as.data.frame(res)

# read in the windows
pass_windows <- readRDS(windowfile)

# read in the G matrix
pass_G <- data.table::fread(grmfile, header = F)
pass_G <- as.matrix(pass_G)

# read in the subset metadata
meta <- as.data.frame(data.table::fread(metafile, header = F))
colnames(meta) <- c("chr", "position")

# clean
gc(); gc();


# run
sims <- sim_gen(x = x, meta = meta, iters = iters, center = T, scheme = "gwas", effect_distribution = rbayesB_fixed,
                parameter_distributions = list(sites = "joint", scale = "joint", d.f = df_func), 
                h_dist = function(x) rnorm(x, hmean, hsd), joint_res = res, joint_acceptance = 0.005,
                peak_delta = .5, peak_pcut = 0.001,
                joint_res_dist = "ks", pass_windows = pass_windows, pass_G = pass_G, GMMAT_infile = gmmatfile)


data.table::fwrite(sims$stats, outname, sep = "\t", quote = F, col.names = F, row.names = F)


