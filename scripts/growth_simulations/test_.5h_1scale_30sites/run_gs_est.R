library(GeneArchEst)
# set parameters
args <- commandArgs(TRUE)
outfile <- as.character(args[1])
hmean <- as.numeric(args[2])
hsd <- as.numeric(args[3])
joint_dist <- readRDS("../../../results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/regression.RDS")
dat <- process_ms("../../../data/theta4k_1000_10_rho40k.txt", chr.length = 10000000, fix_overlaps = T)
genotypes <- dat$x
real.phenotypes <-  unlist(read.table("../../../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.phenos.txt"))
meta <- read.table("../../../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.meta.txt", header = T)

# run an ABC vs the real phenotypes to check fit:
tparam <- sample_joint_quantile(1, joint_dist)
ABC <- ABC_on_hyperparameters(genotypes, phenotypes = real.phenotypes, iters = 1, effect_distribution = rbayesB_fixed,
                              parameter_distributions = list(sites = function(x) rep(tparam$sites, 1),
                                                             d.f = function(x) runif(x, 1, 100),
                                                             scale = function(x) rep(tparam$scale, 1)), 
                              phased = T, save_effects = T, h_dist = function(x) rnorm(x, hmean, hsd), center = T 
)

# get phenotypes
meta$effect <- ABC$effects[1,]
h <- ABC$dist$h
phenotypes <- get.pheno.vals(genotypes, meta$effect, h, phased = T)


# rescale the phenotypes in the first generation so that it has the same shift in optimum over time
scaling_factor <- sd(real.phenotypes)/sd(phenotypes$p)
phenotypes$p <- phenotypes$p*scaling_factor

# real results
est.res <- gs(genotypes = genotypes, meta = meta, phenotypes = phenotypes$p, h = h, gens = 100, var.theta = .1, plot_during_progress = T, chr.length = rep(10000000, 10),
              growth.function = function(n) BL_growth(n, 2),
              survival.function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1.8),
              selection.shift.function = function(opt, iv) optimum_constant_slide(opt, iv, .25),
              K_thin_post_surv = 500) # best immitation of BL's methods.

# concat results
output <- data.frame(method = "est", gen = est.res$run_vars[,"gen"], n = est.res$run_vars[,"N"], ks = ABC$dist$ks) # real

write.table(output, outfile, quote = F, row.names = F, col.names = F)
