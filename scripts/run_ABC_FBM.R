# takes four arguments: input dataset, RDS list containing [[1]] full genotypes, [[2]] genotypes w/ missing,
# [[3]] imputed genotypes, $meta containing snp metadata, and $phenos containing phenotypes;
# name of the output file;
# data_type: specifies if "full" or "imputed" datasets should be used
# h_sd: specifies the sd to use for the h dist


library(GeneArchEst)

args <- commandArgs(TRUE)
dat <- as.character(args[1])
outname <- as.character(args[2])
hsd <- as.numeric(args[3])
hmean <- as.numeric(args[4])

# x <- "../genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS"
# output <- "../genomic_prediction_accuracy/ABC/pi_9999_scale_1_h_5_df_5/ABC_scheme_D_.RDS"

# ABC params
## priors
#pi_func <- function(x) rbeta(x, 20000, 1)
sites_func <- function(n) floor(rexp(n, rate = .0025))
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) rbeta(x, 1, 3)*100
h <- function(x) rnorm(x, hmean, hsd)

## run params
iters <- 1
par <- F

#========================read in genomic data, assign effects, run ABC======
genofile <- paste0(dat, ".rds")
metafile <- paste0(dat, ".geno.meta.RDS")
cat("Genotype RDS:", genofile, "\n")
cat("Metadata and phenotype RDS: ", metafile, "\n")

x <- bigstatsr::big_attach(genofile)
meta <- readRDS(metafile)
phenos <- meta$phenos
meta <- meta$meta
gc(); gc();

# run
ABC_res <- ABC_on_hyperparameters(x = x, phenotypes = phenos$p, iters = iters, 
                                  effect_distribution = rbayesB_fixed, 
                                  h_dist = h,
                                  parameter_distributions = list(sites = sites_func, d.f = df_func, scale = scale_func), 
                                  par = par, center = T)

# save
write.table(ABC_res, outname, col.names = F, row.names = F, quote = F, sep = "\t")

