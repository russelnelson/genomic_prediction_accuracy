# set parameters
args <- commandArgs(TRUE)
outfile <- as.character(args[1])
dat <- process_ms("C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt", chr.length = 10000000, fix_overlaps = T)
genotypes <- dat$x
phenotypes <- unlist(read.table("../genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.phenos.txt"))
meta <- read.table("../genomic_prediction_accuracy/data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.meta.txt", header = T)


# real results
real.res <- gs(genotypes = genotypes, meta = meta, phenotypes = phenotypes, h = .5, gens = 100, var.theta = .1, plot_during_progress = T, chr.length = rep(10000000, 10),
               growth.function = function(n) BL_growth(n, 2),
               survival.function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1.8),
               selection.shift.function = function(opt, iv) optimum_constant_slide(opt, iv, .25),
               K_thin_post_surv = 500) # best immitation of BL's methods.


# BL and breeders
BL.res <- gs_BL(phenotypes, h = .5, K = 500, omega = 1.8, B = 2, var.theta = .1, k = .25, gens = 100)
breeders.res <- gs_breeders(phenotypes, h = .5, K = 500, B = 2, omega = 1.8, var.theta = .1, k = .25, gens = 100)


# genomic prediction
pr_gp <- readRDS("../../results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/pr_gp.RDS")

gp.res <- gs(genotypes[pr_gp$kept.snps,], meta = cbind(meta[pr_gp$kept.snps, 1:2], effect = pr_gp$e.eff$V2),
             phenotypes = phenotypes, h = .5,
             var.theta = .1, gens = 10, plot_during_progress = T, chr.length = rep(10000000, 10),
             growth.function = function(n) BL_growth(n, 2),
             survival.function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1.8),
             selection.shift.function = function(opt, iv) optimum_constant_slide(opt, iv, 0.25),
             K_thin_post_surv = 500)
rm(pr_gp)


# rf
pr_rf <- readRDS("../../results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/pr_rf.RDS")

rf.res <- gs(genotypes[pr_rf$kept.snps,], meta = cbind(meta[pr_rf$kept.snps, 1:2]),
             phenotypes = phenotypes, h = .5,
             var.theta = .1, gens = 100, plot_during_progress = T, chr.length = rep(10000000, 10),
             model = pr_rf$output.model$model,
             growth.function = function(n) BL_growth(n, 2),
             survival.function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1.8),
             selection.shift.function = function(opt, iv) optimum_constant_slide(opt, iv, 0.25),
             K_thin_post_surv = 500)
rm(pr_rf)

# GWAS
pr_gmmat <- pred(genotypes, phenotypes = phenotypes,
                 prediction.program = "GMMAT",
                 maf.filt = F, runID = "gmmat_real")

gmmat.res <- gs_outliers(genotypes, pvals = pr_gmmat$e.eff$PVAL, scores = pr_gmmat$e.eff$SCORE, gens = 100,
                         surivival.function = function(prop_high, opt_prop_high, h.av, ...) prop_survival_scaled_omega(prop_high, opt_prop_high, h.av, omega = 1.8),
                         opt_prop_slide = optimum_proportion_constant_slide(opt, iv, slide = 0.4), pcrit = 1*10^-4)

# concat results
output <- cbind(data.frame(method = "real", gen = real.res$run_vars[,"gen"], n = real.res$run_vars[,"N"]), # real
                data.frame(method = "BL", gen = BL.res$t, n = BL.res$n), # BL
                data.frame(method = "breeders", gen = breeders.res$t, n = breeders.res$n), # breeders
                data.frame(method = "prediction", gen = gp.res$run_vars[,"gen"], n = gp.res$run_vars[,"N"]), # prediction
                data.frame(method = "RF", gen = rf.res$run_vars[,"gen"], n = rf.res$run_vars[,"N"]), # random forest
                data.frame(method = "outliers", gen = gmmat.res$demographics$gen, n = gmmat.res$demographics$n) # gwas outliers
                )

write.table(output, outfile, quote = F, row.names = F, col.names = F)
