library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
h_func <- function(x) rep(.75, x)

args <- commandArgs(TRUE)
x <- as.character(args[1])
y <- as.character(args[2])
outname <- as.character(args[3])
sim_iter <- as.character(args[4])

x <- readRDS("../genomic_prediction_accuracy/ABC/pi_.999_df_5_scale_1_h.75/ABC_input_.999pi_.75h_5df_scale1.RDS")
meta <- x$meta
x <- x$x

res <- readRDS("../genomic_prediction_accuracy/ABC/pi_.999_df_5_scale_1_h.75/ABC_scheme_D_dual_opt_res_pi.999_scale_1.RDS")

sims <- sim_gen(x = x, meta = meta, iters = 1, chr = "group",
                pi_func = "joint", df_func = df_func,
                scale_func = "joint", h_func = h_func, joint_res = res,
                par = F, joint_res_dist = "dist", joint_acceptance = .05)

sims$stats$iter <- sim_iter
write.table(sims$stats, paste0("stats_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(sims$effects, paste0("effects_", sim_iter, "_", outname), sep = "\t", quote = F, col.names = F, row.names = F)


