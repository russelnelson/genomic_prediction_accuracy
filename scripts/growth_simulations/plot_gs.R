library(ggplot2); library(GeneArchEst)

files <- list.files("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/gs/", "gs_est*")
gs_est <- vector("list", length(files))
for(i in 1:length(files)){
  cat(i, "\n")
  gs_est[[i]] <- read.table(paste0("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/gs/", files[i]))
}

gs_est <- dplyr::bind_rows(gs_est)
colnames(gs_est) <- c("method", "generation", "N", "dist", "iter")

best <- 100
gs_sort <- gs_est[,4:5]
gs_sort <- unique(gs_sort)
gs_sort <- gs_sort[order(gs_sort$dist),]
gs_sort <- gs_sort[1:best,]$iter

gs_est_best <- gs_est[which(gs_est$iter %in% gs_sort),]

files <- list.files("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/gs/", "gs_r.+")
gs_other <- vector("list", length(files))
for(i in 1:length(files)){
  cat(i, "\n")
  gs_other[[i]] <- read.table(paste0("results/test_data/sites_30_scale_1_h_.5/hsd.05/trial_1/gs/", files[i]))
}

gs_other <- dplyr::bind_rows(gs_other)
colnames(gs_other) <- c("method", "generation", "N", "iter")
gs_other$dist <- NA

gs_both <- dplyr::bind_rows(list(gs_other, gs_est_best))


gs_both <- gs_both[-which(is.na(gs_both$N)),]






# fix BL res, not stochastic so only need to run once
phenotypes <- unlist(read.table("data/test_data/bayesB_fixed_30_sites_h_.5_scale_1_trial_2.gt.phenos.txt"))
meta <- read.table("data/test_data/bayesB_fixed_30_sites_h_.5_scale_1_trial_2.gt.meta.txt", header = T)
BL.res <- gs_BL(phenotypes, h = .5, K = 250, omega = 1.8, B = 2, var.theta = .1, k = .25, gens = 100)

gs_both <- gs_both[gs_both$method != "BL",]
gs_both <- rbind(gs_both, data.frame(method = "BL", generation = BL.res$t, N = BL.res$n, iter = 1, dist = NA))

ggplot(gs_both, aes(x = generation, y = N, group = iter, color = dist)) + geom_line() +
  facet_wrap(~method) + theme_bw() + scale_color_viridis_c() + theme(strip.background = element_blank())


