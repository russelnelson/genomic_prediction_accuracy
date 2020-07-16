x <- "../genomic_prediction_accuracy/ABC/ABC_input_scale_1_pi_9999_h_5_df_5.RDS"
output <- "../genomic_prediction_accuracy/ABC/pi_9999_scale_1_h_5_df_5/ABC_scheme_D_.RDS"

# ABC params
## priors
pi_func <- function(x) rbeta(x, 200, 1)
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) rbeta(x, 1, 3)*100
h <- function(x) rnorm(x, .5, .1)

## run params
iters <- 100000
par <- 24

#========================read in genomic data, assign effects, run ABC======
d <- readRDS(x)
x <- d[[3]]
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
saveRDS(ABC_res, output)

