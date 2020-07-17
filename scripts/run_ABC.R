# takes four arguments: input dataset, RDS list containing [[1]] full genotypes, [[2]] genotypes w/ missing,
# [[3]] imputed genotypes, $meta containing snp metadata, and $phenos containing phenotypes;
# name of the output file;
# data_type: specifies if "full" or "imputed" datasets should be used
# h_sd: specifies the sd to use for the h dist


library(GeneArchEst)

args <- commandArgs(TRUE)
x <- as.character(args[1])
outname <- as.character(args[2])
data_type <- as.character(args[3])
hsd <- as.numeric(args[4])

# x <- "../genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS"
# output <- "../genomic_prediction_accuracy/ABC/pi_9999_scale_1_h_5_df_5/ABC_scheme_D_.RDS"

# ABC params
## priors
pi_func <- function(x) rbeta(x, 200, 1)
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) rbeta(x, 1, 3)*100
h <- function(x) rnorm(x, .5, hsd)

## run params
iters <- 100000
par <- 24

#========================read in genomic data, assign effects, run ABC======
d <- readRDS(x)
if(data_type == "imputed"){
  x <- d[[3]]
}
if(data_type == "full"){
  x <- d[[1]]
}
meta <- d$meta
phenos <- d$phenos
rm(d); gc(); gc();


# run
ABC_res <- ABC_on_hyperparameters(x = x, phenotypes = phenos$p, iters = iters, 
                                  effect_distribution = rbayesB, 
                                  h_dist = h,
                                  parameter_distributions = list(pi = pi_func, d.f = df_func, scale = scale_func), 
                                  par = par, center = T)

# save
saveRDS(ABC_res, outname)

