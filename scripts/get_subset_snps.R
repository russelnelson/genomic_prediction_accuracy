library(bigstatsr)
# import data
dat <- big_attach(rdsfile = "/home/hemstrow/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.rds")

# find mafs
rs <- big_apply(dat, a.FUN = function(X, ind) rowSums(X[, ind]),a.combine ='plus')
maf <- rs/(2*ncol(dat))
too.big <- which(maf > 0.5)
maf[too.big] <- 1 - maf[too.big]

# get snps to keep
good <- which(maf >= 0.05)
sub_snps <- sort(sample(good, 1000000))

# save
write.table(sub_snps, "/home/hemstrow/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/subset_snps.txt", row.names = F, col.names = F, quote = F)

