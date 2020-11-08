dat <- process_ms("C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt", chr.length = 10000000, fix_overlaps = T)

meta <- dat$meta
meta$effect <- rbayesB_fixed(nrow(dat$x), sites = 30, d.f = 5, scale = 1)

# assign phenotypes
phenos <- get.pheno.vals(dat$x, meta$effect, .5, phased = T)

# add missing data
x_missing <- sprinkle_missing_data(dat$x)
x_imp <- impute_and_phase_beagle(x_missing, meta, "D://usr/bin/beagle.18May20.d20.jar", num_threads = 8)

write.table(meta, "../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.meta.txt", quote = F, row.names = F, col.names = T)
write.table(data.frame(p = phenos$p), "../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.phenos.txt", quote = F, row.names = F, col.names = F)
file.copy("data.gt.vcf.gz", "../data/test_data/bayesB_fixed_30_sites_h_.5_scale_1.gt.vcf")
