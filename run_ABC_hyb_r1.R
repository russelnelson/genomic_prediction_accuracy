x <- "ABC/ABC_input_.9999pi_.75h_5df_scale1_r2.RDS"
output <- "ABC/pi_.9999_df_5_scale_1_h.75/ABC_scheme_D_dual_opt_res_pi.9999_scale_1_dist_trial_r2.RDS"

# ABC params
## priors
pi_func <- function(x) rbeta(x, 200, 1)
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) rbeta(x, 1, 3)*100
h <- 0.75

## run params
iters <- 100000
scheme <- "D"
par <- 6

#========================source scripts===============================
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.4", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(methods); library(data.table); library(doParallel); library(dplyr); library(GenKern)

source("~/coalescence/prediction_accuracy/genomic_prediction_accuracy/growth_sim.R")

#========================read in genomic data, assign effects, run ABC======

d <- readRDS(x)
x <- d$x
meta <- d$meta
phenos <- d$phenos
rm(d); gc(); gc();

# run
ABC_res <- ABC_on_hyperparameters(x, phenos$p, iters = iters, pi_func = pi_func, 
                                  df_func = df_func, scale_func = scale_func, ABC_scheme = "D",
                                  h = h, par = par, save_effects = F)

# save
saveRDS(ABC_res$dists, output)

