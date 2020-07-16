x <- "C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt"
outname <- "../genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS"

# read in data
x <- process_ms(x, 10000000)
meta <- x$meta
x <- x$x

# assign effects
meta$effect <- rbayesB(nrow(meta), pi = .9999, d.f = 5, scale = 1)

# assign phenotypes
phenos <- get.pheno.vals(x, meta$effect, .5)

# add missing data
x_missing <- sprinkle_missing_data(x)
x_imp <- impute_and_phase_beagle(x_missing, meta, "D://usr/bin/beagle.18May20.d20.jar", num_threads = 8)
# can assess imputation success with assess_imputation if wanted.

saveRDS(list(x = x, x_missing =  x_missing, x_imp = x_imp, meta = meta, phenos = phenos), outname)
