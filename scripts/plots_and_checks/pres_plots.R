library(GeneArchEst); library(ggplot2)
#================30 and 300 sites test runs==============
# 30 sites

rf30 <- readRDS("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/RF.RDS")
reg30 <- readRDS("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/regression.RDS")
abc30 <- data.table::fread("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/ABC_res.txt")
colnames(abc30) <- c("sites", "d.f", "scale", "h", GeneArchEst::names_diff_stats)
abc30$hits <- ifelse(abc30$ks <= quantile(abc30$ks, 0.005, na.rm = T), 1, 0)

rf30p <- plot_rf(rf30, "sites")
reg30p <- plot_reg(reg30)
reg30p <- reg30p$joint_quantile_plot + geom_vline(aes(xintercept = 30)) + 
  geom_hline(aes(yintercept = 1))
abc30p <- ggplot(abc30[abc30$hits == 1,], aes(x = sites, y = scale, color = ks)) +
  geom_point() + theme_bw() +  
  geom_point(data = abc30[abc30$hits == 0,], alpha = .01, color = "grey") +
  scale_color_viridis_c()

# 300 sites
rf300 <- readRDS("results/test_data/sites_300_scale_1_h_.5/hsd.05/RF.RDS")
reg300 <- readRDS("results/test_data/sites_300_scale_1_h_.5/hsd.05/regression.RDS")
abc300 <- data.table::fread("results/test_data/sites_300_scale_1_h_.5/hsd.05/ABC_res.txt")
colnames(abc300) <- c("sites", "d.f", "scale", "h", GeneArchEst::names_diff_stats)
abc300$hits <- ifelse(abc300$ks <= quantile(abc300$ks, 0.005, na.rm = T), 1, 0)

rf300p <- plot_rf(rf300, "sites")
reg300p <- plot_reg(reg300)
reg300p <- reg300p$joint_quantile_plot + geom_vline(aes(xintercept = 300)) + 
  geom_hline(aes(yintercept = 1))
abc300p <- ggplot(abc300[abc300$hits == 1,], aes(x = sites, y = scale, color = ks)) +
  geom_point() + theme_bw() +  
  geom_point(data = abc300[abc300$hits == 0,], alpha = .01, color = "grey") +
  scale_color_viridis_c()


gridExtra::grid.arrange(reg30p + theme(legend.position = "none"), rf30p$cross_val + theme(legend.position = "none"), rf30p$density,
                        layout_matrix = matrix(c(1, 2, 1, 3), nrow = 2))

gridExtra::grid.arrange(reg300p + theme(legend.position = "none"), rf300p$cross_val + theme(legend.position = "none"), rf300p$density,
                        layout_matrix = matrix(c(1, 2, 1, 3), nrow = 2))

gridExtra::grid.arrange(reg30p, reg300p, rf30p$cross_val, rf300p$cross_val,
                        nrow = 2, ncol = 2)