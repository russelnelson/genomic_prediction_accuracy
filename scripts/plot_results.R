library(GeneArchEst); library(ggplot2)
rf <- readRDS("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/RF.RDS")
reg <- readRDS("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/regression.RDS")
abc <- read.table("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/ABC_res.txt")
sims <- read.table("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/sim_gen_stats.txt")

colnames(abc)[1:2] <- c("sites", "scale")
colnames(sims)[1:2] <- c("sites", "scale")


rfp <- plot_rf(rf, hyperparameter = "sites")
regp <- plot_reg(reg)
regp$joint_quantile_plot +
  #geom_point(data = sims, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.05) +
  geom_point(data = abc, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.5) + scale_y_continuous(limits = c(0, 5))
  

rf$sites$point_estimate
reg$optimal_fits

lm <- cbind(c(1, 1, 3), c(2, 2, 3))
grid.arrange(regp$joint_quantile_plot, rfp$density, rfp$cross_val, ncol = 2, layout_matrix = lm)
