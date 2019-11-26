args <- commandArgs(TRUE)
x <- as.character(args[1])
outname <- as.character(args[2])
run_n <- as.character(args[3])

# file and storage information:
runID <- paste0("R_", run_n) # run ID and directory name where intermediate results will be stored. Will be created if needed.
save.meta <- T # should we save the metadata file when doing prediction/GWAS

# genome information
chrl <- 10000000 # chromosome length CURRENTLY MUST BE THE SAME ACROSS ALL CHRs

# effect size information
effect.dist <- "bayesB" # which effect distribution should we use
prob.effect <- 0.01 # probability that any single SNP has an effect, for effect dists without a fixed number of effect loci.
h <- 0.75 # h^2, or heritability for the trait prior to any selection.
pi <- 1 - prob.effect # probability that a site has zero effect for "bayesB"
d.f <- 5 # degrees of freedom for the "bayesB" t-distribution
scale <- 1 # scale factor for the "bayesB" t-distribution

# chain information
chain_length <- 100000
burnin <- 5000
thin <- 100

# ABC params
## priors
pi_func <- function(x) rbeta(x, 25, 1)
df_func <- function(x) runif(x, 3, 7)
scale_func <- function(x) runif(x, 1, 1)

#========================source scripts===============================
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.4", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(methods); library(data.table)

source("~/coalescence/prediction_accuracy/genomic_prediction_accuracy/growth_sim.R")

#========================read in genomic data, assign effects, run ABC======



# read in data
x <- process_ms(x, chrl)
meta <- x$meta
x <- x$x

# assign effects
meta$effect <- rbayesB(nrow(meta), pi, d.f, scale)
cat("Mean effect size:", mean(meta$effect), "\n")

# assign phenotypes
phenos <- get.pheno.vals(x, meta$effect, .75)

# run
ABC_res <- ABC_on_hyperparameters(x, phenos$p, iters = 1, pi_func = pi_func, 
                               df_func = df_func, scale_func = scale_func, ABC_scheme = "C",
                               h = h, chain_length = chain_length, burnin = burnin, thin = thin, run_number = run_n)

# save
write.table(ABC_res, paste0("./", outname, "_", runID, "_dist.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
