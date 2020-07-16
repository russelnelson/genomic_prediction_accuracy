library(snpR); library(ggplot2)

#==========40kr=========
infile <- "theta4k_1000_1_rho40k.txt" #infile
LD_chr1 <- LD_full_pairwise(infile, NULL, T, F, input = "ms", chr.length = 10000000, par = 3, maf = 0.05)
LD_chr1$p1 <- round(LD_chr1$p1)
LD_chr1$p2 <- round(LD_chr1$p2)
LD_chr1$proximity <- abs(LD_chr1$p1 - LD_chr1$p2)
saveRDS(LD_chr1, "LD.out.r40k.chr1.RDS")


ss <- sample(nrow(LD_chr1), 10000, F)
ggplot(LD_chr1[ss,], aes(x = proximity, y = rsq)) + geom_point() + theme_bw() + geom_smooth()
ggplot(LD_chr1[ss,], aes(x = proximity, y = Dprime)) + geom_point() + theme_bw() + geom_smooth()

ggplot(LD_chr1[LD_chr1$proximity <= 50000,], aes(x = proximity, y = rsq)) + 
  geom_point(alpha = 0.05) + 
  theme_bw() + geom_smooth()
ggplot(LD_chr1[LD_chr1$proximity <= 50000,], aes(x = proximity, y = Dprime)) + 
  geom_point(alpha = 0.05) +
  theme_bw() + geom_smooth()

#==========8kr
LD_chr1 <- readRDS("LD.out.ss.chr1.RDS")
ggplot(LD_chr1[LD_chr1$proximity <= 50000,], aes(x = proximity, y = rsq)) + 
  geom_point(alpha = 0.05) + 
  theme_bw() + geom_smooth()
ggplot(LD_chr1[LD_chr1$proximity <= 50000,], aes(x = proximity, y = Dprime)) + 
  geom_point(alpha = 0.05) +
  theme_bw() + geom_smooth()