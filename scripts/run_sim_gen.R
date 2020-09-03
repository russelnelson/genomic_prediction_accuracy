library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
df_func <- function(x) runif(x, 1, 100)


args <- commandArgs(TRUE)
x <- as.character(args[1])
ABCfile <- as.character(args[2])
outname <- as.character(args[3])
data_type <- as.character(args[4])

iters <- 1

x <- readRDS(x)
meta <- x$meta
if(data_type == "imputed"){
  x <- x[[3]]
}
if(data_type == "full"){
  x <- x[[1]]
}

res <- data.table::fread(ABCfile, header = F)
colnames(res) <- c("pi", "d.f", "scale", "h", GeneArchEst::names_diff_stats)
res <- as.data.frame(res)

sims <- sim_gen(x = x, meta = meta, iters = iters, center = T, scheme = "gwas", 
                parameter_distributions = list(pi = "joint", scale = "joint", d.f = df_func), 
                h_dist = function(x) rnorm(x, .5, .1), joint_res = res, joint_acceptance = 0.005, 
                joint_res_dist = "ks")


data.table::fwrite(sims$stats, outname, sep = "\t", quote = F, col.names = F, row.names = F)


