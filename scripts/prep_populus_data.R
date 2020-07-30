# read and prep the populus data
## read
dat <- process_vcf("~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/gatk_882_WG_genotypes_biallelic_snps_VQSR_non_filtered_bis.gt.vcf.gz")
phenos <- data.table::fread("~/coalescence/prediction_accuracy/genomic_prediction_accuracy/data/populus_data/sorted_phenotypes.txt")
meta <- dat$meta

## clean and save
dat <- clean_phenotypes(dat$genotypes, phenos$`Stomatal density (SD)`)
saveRDS(list(full = FALSE, missing = FALSE, imputed = dat$genotypes, meta = meta, phenos = list(p = dat$phenos)),
        "~/coalescence/prediction_accuracy/genomic_prediction_accuracy/results/populus/geno_pheno.RDS")
