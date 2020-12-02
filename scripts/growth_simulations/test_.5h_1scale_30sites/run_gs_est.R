library(GeneArchEst)
# set parameters
args <- commandArgs(TRUE)
outfile <- as.character(args[1])
dat <- process_ms("C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt", chr.length = 10000000, fix_overlaps = T)
genotypes <- dat$x
meta <- read.table("../genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.meta.txt", header = T)
joint_dist <- readRDS("../genomic_prediction_accuracy/results/test_data/sites_30_scale_1_h_.5/hsd.05/regression.RDS")
real.phenotypes <- unlist(read.table("../../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.phenos.txt"))

tparam <- sample_joint_quantile(1, joint_dist)
meta$effect <- rbayesB_fixed(nrow(meta), sites = tparam$sites, scale = tparam$scale, d.f = 5)
h <- rnorm(1, .5, .05)
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
output <- data.frame(method = "est", gen = est.res$run_vars[,"gen"], n = est.res$run_vars[,"N"]), # real
                
write.table(output, outfile, quote = F, row.names = F, col.names = F)
