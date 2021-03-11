# takes three arguments: input dataset, RDS list containing $genos full genotypes, $meta containing snp metadata, and $phenos containing phenotypes;
# name of the output file;
# data_type: specifies if "full" or "imputed" datasets should be used
# h_sd: specifies the sd to use for the h dist


library(GeneArchEst)

args <- commandArgs(TRUE)
x <- as.character(args[1])
outname <- as.character(args[2])
hsd <- as.numeric(args[3])
hmean <- as.numeric(args[4])

# x <- "../genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS"
# output <- "../genomic_prediction_accuracy/ABC/pi_9999_scale_1_h_5_df_5/ABC_scheme_D_.RDS"

# ABC params
## priors
sites_func <- function(n) floor(runif(n, 0, 1000))
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) runif(x, 0, 100)
h <- function(x) rnorm(x, hmean, hsd)



## run params
iters <- 100000
par <- 24

#========================read in genomic data, assign effects, run ABC======

x <- readRDS(x)
phenos <- x$phenos
meta <- x$meta
meta <- as.data.frame(meta)
x <- x$genos

# run
ABC_res <- ABC_on_hyperparameters(x = x, phenotypes = phenos$p, iters = iters, 
                                  effect_distribution = rbayesB_fixed, 
                                  h_dist = h,
                                  parameter_distributions = list(sites = sites_func, d.f = df_func, scale = scale_func), 
                                  par = par, center = T, phased = F)

# save
#saveRDS(ABC_res, outname)
write.table(ABC_res$ABC_res, outname, col.names = F, row.names = F, quote = F, sep = "\t")


