library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
df_func <- function(x) runif(x, 1, 100)


args <- commandArgs(TRUE)
x <- as.character(args[1])
ABCfile <- as.character(args[2])
outname <- as.character(args[3])
hsd <- as.numeric(args[4])
hmean <- as.numeric(args[5])
peak_delta <- as.numeric(args[6])
peak_pcut <- as.numeric(args[7])

iters <- 1

x <- readRDS(x)
phenos <- x$phenos
meta <- x$meta
meta <- as.data.frame(meta)
x <- x$genos

res <- data.table::fread(ABCfile, header = F)
colnames(res) <- c("sites", "d.f", "scale", "h", GeneArchEst::names_diff_stats)
res <- as.data.frame(res)


sims <- sim_gen(x = x, meta = meta, iters = iters, center = T, scheme = "gwas", effect_distribution = rbayesB_fixed, 
                parameter_distributions = list(sites = "joint", scale = "joint", d.f = df_func), 
                h_dist = function(x) rnorm(x, hmean, hsd), joint_res = res, joint_acceptance = 0.005,
                peak_delta = peak_delta, peak_pcut = peak_pcut, 
                joint_res_dist = "ks", phased = FALSE)

data.table::fwrite(sims$stats, outname, sep = "\t", quote = F, col.names = F, row.names = F)


