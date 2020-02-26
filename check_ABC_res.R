# import data
res <- read.table("ABC/cat_dists_ABC_.9999pi_.75h_5df_scale1_scale_not_fixed_hyb_out_gmmat2")

#colnames(res) <- c("pi", "df", "scale", "dist", "norm.dist", "iter")
comp_cols <- c("kurtosis", "skewness", "mean", "median", "sd", "var",
               paste0("Quantile_", seq(0.1, 0.9, length.out = 20)),
               "ks", "norm.ks", "lepage", "cucconi")

colnames(res) <- c("pi", "df", "scale", comp_cols, "peak_count_diff", paste0("peak_", comp_cols), "iter")

res <- as.data.frame(res)


# find the optimum pi and scale via gam smoothing. The difference in the number of peaks
# seems to relate only to pi, so do that first, then use the fit optimum pi to predict scale, since the two
# are tightly linked. For the second pass, use only the runs with the best ks scores.
library(ggplot2); library(mgcv)
# pi gam
p <- ggplot(res, aes(x = log10(1 - pi), y = peak_count_diff)) + theme_bw() +
  geom_smooth() + geom_point(alpha = .25)
res$trans_pi <- log10(1 - res$pi)
g <- gam(peak_count_diff ~ s(trans_pi, bs = "cs", k = 100), data = res)

# pred pi gam
pd <- data.frame(pi = seq(min(res$pi), max(res$pi), by = .000001))
pd$trans_pi <- log10(1 - pd$pi)
pd$pred_pcd <- unlist(predict(g, newdata = pd))
opt.pi <- pd[which.min(pd$pred_pcd),]$pi


# scale gam
ores <- res
res <- res[res$ks <= quantile(res$ks, .01),]
g2 <- gam(scale ~ s(trans_pi, bs = "cs"), data = res)

# pred scale gam
pd2 <- data.frame(trans_pi = log10(1 - opt.pi))
opt.scale <- unlist(predict(g2, newdata = pd2))







# res$hits <- ifelse(res$lepage <= quantile(res$lepage, .01), 1, 0)
# res$hits <- ifelse(res$peak_count_diff <= 1, 1, 0)
# ores <- res
# res <- res[res$hits == 1,]
# res$hits <- ifelse(res$ks <= quantile(res$ks, .1) &
#                      res$peak_ks <= quantile(res$peak_ks, .1), 1, 0)

# res$norm.hits <- ifelse(res$norm.dist <= quantile(res$norm.dist, .01, na.rm = T), 1, 0)
# pull the best pi value using a gam smooth


# 
# ggplot(res, aes(x = pi, y = df, color = as.factor(hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# ggplot(res, aes(x = pi, y = ks, color = as.factor(hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# 
# ggplot(res, aes(x = df, y = ks, color = as.factor(hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# ggplot(res, aes(x = scale, y = ks, color = as.factor(hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# ggplot(res, aes(x = pi, y = scale, color = as.factor(hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# 
# # 
# # ggplot(res, aes(x = pi, y = df, color = as.factor(norm.hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# # ggplot(res, aes(x = pi, y = norm.dist, color = as.factor(norm.hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# # 
# # ggplot(res, aes(x = df, y = norm.dist, color = as.factor(norm.hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# # ggplot(res, aes(x = scale, y = norm.dist, color = as.factor(norm.hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# # ggplot(res, aes(x = scale, y = pi, color = as.factor(norm.hits))) + geom_point() + scale_color_viridis_d() + theme_bw()
# # 
# ggplot(res, aes(x = pi, y = ks, color = as.factor(hits))) + geom_point() + 
#   scale_color_viridis_d() + theme_bw() + scale_y_continuous(limits = c(0, 0.005))
# 
# ggplot(res, aes(x = pi, y = scale, z = hits)) + stat_summary_hex(bins = 20) + scale_fill_viridis_c() + theme_bw()
# ggplot(res, aes(x = pi, y = scale, z = peak_count_diff)) + stat_summary_hex(bins = 20) + scale_fill_viridis_c() + theme_bw()
# 
# ggplot(res[res$hits == 1,], aes(x = pi, y = scale, color = ks)) + geom_point() + scale_color_viridis_c()
# 
# gres <- res[res$hits == 1,]
# mean(gres$pi)
# mean(gres$df)
# mean(gres$scale)
# 
# 
# library(readr)
# p.gs.res <- read_delim("ABC/cat_gs_pseudo_ABC_.9999pi_.75h_5df_scale1_scale_not_fixed_hyb_out.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
# colnames(p.gs.res) <- c("run", "N", "mu_pheno", "mu_a", "opt", "diff", "var_a", "stochastic_opt", "gen", "source", "iter")
# r.gs.res <- read_delim("ABC/cat_gs_real_ABC_.9999pi_.75h_5df_scale1_scale_not_fixed_hyb_out.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
# colnames(r.gs.res) <- colnames(p.gs.res)
# r.gs.res$rID <- paste0(r.gs.res$iter, "_", r.gs.res$run)
# 
# p2 <- plot_sim_res(res, p.gs.res, real.sims = r.gs.res[r.gs.res$iter < 100,], acceptance_ratio = .01,
#                    smooth = T)
# p <- plot_sim_res(res, p.gs.res, real.sims = r.gs.res[r.gs.res$iter < 100,], acceptance_ratio = 1, smooth = F)
