dat <- process_ms("../../data/theta4k_1000_10_rho40k.txt", chr.length = 10000000, fix_overlaps = T)
genotypes <- dat$x
phenotypes <- unlist(read.table("../../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.phenos.txt"))
meta <- read.table("../../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.meta.txt", header = T)


# prediction
pr_gp <- pred(genotypes, meta, phenotypes = phenotypes, prediction.program = "BGLR", prediction.model = "BayesB", chain_length = 10000, burnin = 1000,  thin = 100, maf.filt = F, par = 3)

pr_rf <- pred(genotypes, meta, phenotypes = phenotypes, prediction.program = "ranger", ntree = 500, prediction.model = "RJ", maf.filt = F, par = 3)

pr_gmmat <- pred(genotypes, phenotypes = phenotypes,
                 prediction.program = "GMMAT",
                 maf.filt = F, runID = "gmmat_real", par = 3)

# save
saveRDS(pr_gp, "../../results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/pr_gp.RDS")
saveRDS(pr_rf, "../../results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/pr_rf.RDS")
saveRDS(pr_gmmat, "../../results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/pr_gmmat.RDS")

