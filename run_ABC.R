args <- commandArgs(TRUE)
x <- as.character(args[1])
#x <- "C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt"
outname <- as.character(args[2])

# file and storage information:
save.meta <- F # should we save the metadata file when doing prediction/GWAS

# genome information
chrl <- 10000000 # chromosome length CURRENTLY MUST BE THE SAME ACROSS ALL CHRs

# effect size information
effect.dist <- "bayesB" # which effect distribution should we use
prob.effect <- 0.0001 # probability that any single SNP has an effect, for effect dists without a fixed number of effect loci.
h <- 0.75 # h^2, or heritability for the trait prior to any selection.
pi <- 1 - prob.effect # probability that a site has zero effect for "bayesB"
d.f <- 5 # degrees of freedom for the "bayesB" t-distribution
scale <- 1 # scale factor for the "bayesB" t-distribution

# chain information
chain_length <- 10000
burnin <- 500
thin <- 100

# ABC params
## priors
pi_func <- function(x) rbeta(x, 200, 1)
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) runif(x, 1, 1)

#========================source scripts===============================
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.4", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(methods); library(data.table); library(doParallel)

source("~/coalescence/prediction_accuracy/genomic_prediction_accuracy/growth_sim.R")

#========================read in genomic data, assign effects, run ABC======



# read in data
x <- process_ms(x, chrl)
meta <- x$meta
x <- x$x

# assign effects
meta$effect <- rbayesB(nrow(meta), pi, d.f, scale)
cat("Mean effect size:", mean(meta$effect), "\n")

# assign phenotypes
phenos <- get.pheno.vals(x, meta$effect, .75)

# run
ABC_res_pi <- ABC_on_hyperparameters(x, phenos$p, iters = 1, pi_func = pi_func, 
                                     df_func = df_func, scale_func = scale_func, ABC_scheme = "C",
                                     h = h, chain_length = chain_length, burnin = burnin, thin = thin, par = F)

# res <- ABC_res_pi_.9999$dists
# 
# res$hits <- ifelse(res$dist <= quantile(res$dist, .0025), 1, 0)
# gres <- res[res$hits == 1,]
# pseudo.meta <- meta
# best.res <- which.min(res$dist)
# pseudo.meta$effect <- ABC_res_pi_.9999$effects[[best.res]]
# pseudo.phenos <- get.pheno.vals(x, pseudo.meta$effect, h = .75)
# 
# # real
# sim.real.x <- list(e.eff = data.frame(site = 1:nrow(meta), effect = meta$effect), h = h, 
#                    phenotypes = phenos, meta = meta, x = x)
# 
# sim.real <- gs(x = sim.real.x, 
#                gens = 100,
#                pred.method = "real",
#                growth.function =  l_g_func, 
#                survival.function = survival.dist.func, 
#                selection.shift.function = sopt.func, 
#                rec.dist = rec.dist, 
#                var.theta = var.theta,
#                plot_during_progress = F, 
#                chr.length = chrl, 
#                print.all.freqs = F,
#                fgen.pheno = T)
# 
# 
# # pseudo
# sim.psuedo.x <- list(e.eff = data.frame(site = 1:nrow(meta), effect = pseudo.meta$effect), h = h, 
#                      phenotypes = pseudo.phenos, meta = pseudo.meta, x = x)
# 
# sim.pseudo <- gs(x = sim.psuedo.x, 
#                gens = 100,
#                pred.method = "real",
#                growth.function =  l_g_func, 
#                survival.function = survival.dist.func, 
#                selection.shift.function = sopt.func, 
#                rec.dist = rec.dist, 
#                var.theta = var.theta,
#                plot_during_progress = F, 
#                chr.length = chrl, 
#                print.all.freqs = F,
#                fgen.pheno = T)
# 
# 
# 
# comb <- rbind(cbind(as.data.frame(sim.pseudo$run_vars, stringsAsFactors = F), model = "pred", gen = 1:nrow(sim.pseudo$run_vars)),
#               cbind(as.data.frame(sim.real$run_vars, stringsAsFactors = F), model = "real", gen = 1:nrow(sim.real$run_vars)))
# comb$model <- as.character(comb$model)

# save
write.table(ABC_res$dists, paste0("dists_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ABC_res$effects, paste0("effects_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
