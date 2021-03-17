library(GeneArchEst); library(ggplot2)
rf <- readRDS("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/RF.RDS")
reg <- readRDS("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/regression.RDS")
abc <- read.table("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/ABC_res.txt")
sims <- read.table("results/populus/h_0.05_sd_rbayesb_fixed/stomatal_density/sim_gen_stats.txt.gz")

colnames(abc) <- c("sites", "df", "scale", "h", GeneArchEst::names_diff_stats)
colnames(sims)[1:3] <- c("sites", "df", "scale")


rfp <- plot_rf(rf, hyperparameter = "sites")
regp <- plot_reg(reg)
regp$joint_quantile_plot +
  #geom_point(data = sims, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.05) +
  geom_point(data = abc, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.5) + scale_y_continuous(limits = c(0, 5))
  

rf$sites$point_estimate
reg$optimal_fits

lm <- cbind(c(1, 1, 3), c(2, 2, 3))
grid.arrange(regp$joint_quantile_plot, rfp$density, rfp$cross_val, ncol = 2, layout_matrix = lm)

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = sites, y = scale, color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c() +
  geom_point(data = abc, alpha = 0.05, color = "grey")

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = sites, y = scale, color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c() #+
  #geom_point(data = abc, alpha = 0.05, color = "grey")

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = log10(sites), y = log10(scale), color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c()






###################
library(GeneArchEst); library(ggplot2); library(gridExtra)
rf <- readRDS("results/populus/h_0.05_sd_rbayesb_fixed/petiole_length/RF.RDS")
reg <- readRDS("results/populus/h_0.05_sd_rbayesb_fixed/petiole_length/regression.RDS")
abc <- read.table("results/populus/h_0.05_sd_rbayesb_fixed/petiole_length/ABC_res.txt")
sims <- read.table("results/populus/h_0.05_sd_rbayesb_fixed/petiole_length/sim_gen_stats.txt.gz")

colnames(abc) <- c("sites", "df", "scale", "h", GeneArchEst::names_diff_stats)
colnames(sims)[1:3] <- c("sites", "df", "scale")


rfp <- plot_rf(rf, hyperparameter = "sites")
regp <- plot_reg(reg)
regp$joint_quantile_plot +
  #geom_point(data = sims, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.05) +
  geom_point(data = abc, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.5) + scale_y_continuous(limits = c(0, 5))


rf$sites$point_estimate
reg$optimal_fits

lm <- cbind(c(1, 1, 3), c(2, 2, 3))
grid.arrange(regp$joint_quantile_plot, rfp$density, rfp$cross_val, ncol = 2, layout_matrix = lm)

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = sites, y = scale, color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c() +
  geom_point(data = abc, alpha = 0.05, color = "grey")

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = sites, y = scale, color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c() #+
#geom_point(data = abc, alpha = 0.05, color = "grey")

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = log10(sites), y = log10(scale), color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c()



###################
library(GeneArchEst); library(ggplot2); library(gridExtra)
rf <- readRDS("results/test_data/sites_30_scale_1_h_.5/hsd.05/RF.RDS")
reg <- readRDS("results/test_data/sites_30_scale_1_h_.5/hsd.05/regression.RDS")
abc <- read.table("results/test_data/sites_30_scale_1_h_.5/hsd.05/ABC_res.txt")
sims <- read.table("results/test_data/sites_30_scale_1_h_.5/hsd.05/sim_gen_stats.txt.gz")

colnames(abc) <- c("sites", "df", "scale", "h", GeneArchEst::names_diff_stats)
colnames(sims)[1:3] <- c("sites", "df", "scale")


rfp <- plot_rf(rf, hyperparameter = "sites")
regp <- plot_reg(reg)
regp$joint_quantile_plot +
  #geom_point(data = sims, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.05) +
  geom_point(data = abc, mapping = aes(x = sites, y = scale, z = NULL), color = "grey", alpha = 0.5) + scale_y_continuous(limits = c(0, 5))


rf$sites$point_estimate
reg$optimal_fits

lm <- cbind(c(1, 1, 3), c(2, 2, 3))
grid.arrange(regp$joint_quantile_plot, rfp$density, rfp$cross_val, ncol = 2, layout_matrix = lm)

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = sites, y = scale, color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c() +
  geom_point(data = abc, alpha = 0.05, color = "grey")

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = sites, y = scale, color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c() #+
#geom_point(data = abc, alpha = 0.05, color = "grey")

ggplot(abc[which(abc$ks <= quantile(abc$ks, 0.005, na.rm = T)),], aes(x = log10(sites), y = log10(scale), color = ks)) + geom_point() +
  theme_bw() + scale_color_viridis_c()




test <- function(n) floor(rexp(n, rate = .01))
plot(density(test(1000)))


