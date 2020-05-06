library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
h_func <- function(x) rep(.75, x)
df_func <- function(x) runif(x, 1, 100)


args <- commandArgs(TRUE)
x <- as.character(args[1])
y <- as.character(args[2])
outname <- as.character(args[3])
sim_iter <- as.character(args[4])

x <- readRDS(x)
meta <- x$meta
x <- x$x

res <- readRDS(y)

sims <- sim_gen(x = x, meta = meta, iters = 1, chr = "group", center = T,
                pi_func = "joint", df_func = df_func,
                scale_func = "joint", h_func = h_func, joint_res = res,
                par = F, joint_res_dist = "dist", joint_acceptance = .05)

sims$stats$iter <- sim_iter
write.table(sims$stats, paste0("stats_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(sims$effects, paste0("effects_", sim_iter, "_", outname), sep = "\t", quote = F, col.names = F, row.names = F)


