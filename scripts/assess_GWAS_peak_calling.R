phenotypes <- readRDS("../genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.geno.meta.RDS")
phenotypes <- phenotypes$phenos$p
G <- read.table("../genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.subset.Gmat")
G <- as.matrix(G)


gres <- pred_gwas_FBM(phenotypes = phenotypes, pass_G = G,
                      GMMAT_infile = "../genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.subset.012.clean.transpose")


meta <- read.table("../genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.subset.meta")
colnames(meta) <- c("chr", "position")
meta$effect <- gres$e.eff$PVAL
meta$snpID <- 1:nrow(meta)

peaks <- findpeaks_multi(x = meta, delta = .5, pcut = .001, pvals = T)
peaks$chr_pos <- paste0(peaks$chr, "_", peaks$pos)
meta$chr_pos <- paste0(meta$chr, "_", meta$position)
anno.ids <- meta$snpID[which(meta$chr_pos %in% peaks$chr_pos)]

library(snpR)
p <- plot_manhattan(meta, plot_var = "effect", bp = "position", highlight = anno.ids, log.p = T, )
p$plot
