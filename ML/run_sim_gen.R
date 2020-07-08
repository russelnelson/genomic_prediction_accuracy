library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
df_func <- function(x) runif(x, 1, 100)

args <- commandArgs(TRUE)
x <- as.character(args[1])
y <- as.character(args[2])
outname <- as.character(args[3])
par <- 24
iters <- 100000

x <- readRDS(x)
meta <- x$meta
x <- x[[3]]

res <- readRDS(y)$ABC_res

sims <- sim_gen(x = x, meta = meta, iters = iters, center = T, scheme = "gwas", 
                parameter_distributions = list(pi = "joint", scale = "joint", d.f = df_func), 
                h_dist = function(x) rnorm(x, .5, .1), par = par, joint_res = res, joint_acceptance = 0.005, 
                joint_res_dist = "ks")

data.table::fwrite(sims$stats, paste0("stats_", outname), sep = "\t", quote = F, col.names = T, row.names = F)


