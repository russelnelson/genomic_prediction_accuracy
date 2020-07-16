sim_dat <- readr::read_delim("D:/Users/hemst/Documents/GitHub/genomic_prediction_accuracy/ML/sim_stats_scale_1_pi_.9999_h_5_hdist_imp.txt",
                             "\t", escape_double = FALSE, col_names = FALSE,
                             trim_ws = TRUE)
colnames(sim_dat) <- c("pi", "scale", "d.f", "h", names_descriptive_stats)
x <- readRDS("../genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS")
phenos <- x$phenos$p
meta <- x$meta
x <- x[[3]]
res <- readRDS("../genomic_prediction_accuracy/ABC/pi_9999_scale_1_h_5_df_5/ABC_scheme_D_.RDS")

# note, not naming and removing the iter column, and not checking for that error, produces an uninformative error
# in the quantforresterror. Might gather some example data and post it as an issue.
pred <- hyperparameter_random_forest(x = x, meta = meta, phenos = phenos, num_trees = 10000,
                                    sims = sim_dat, num_threads = 4, hyperparameter_to_estimate = "pi")


saveRDS(pred, "../genomic_prediction_accuracy/ML/pi_9999_h_5_df_5_scale_1_hdist_imp_RF.RDS")
