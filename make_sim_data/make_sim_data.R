x <- "C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt"
outname <- "ABC/ABC_input_.999pi_.75h_5df_scale1.RDS"

# genome information
chrl <- 10000000 # chromosome length CURRENTLY MUST BE THE SAME ACROSS ALL CHRs

# effect size information
effect.dist <- "bayesB" # which effect distribution should we use
prob.effect <- 0.001 # probability that any single SNP has an effect, for effect dists without a fixed number of effect loci.
h <- 0.75 # h^2, or heritability for the trait prior to any selection.
pi <- 1 - prob.effect # probability that a site has zero effect for "bayesB"
d.f <- 5 # degrees of freedom for the "bayesB" t-distribution
scale <- 1 # scale factor for the "bayesB" t-distribution

#========================source scripts===============================
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.4", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(methods); library(data.table); library(doParallel); library(dplyr)

source("growth_sim.R")


#========================models and distributions for simulations=======================================
effect.dist.func.znorm <- function(n){
  eff <- rbinom(n, 1, prob.effect) #does each site have an effect?
  eff[which(eff == 1)] <- rnorm(sum(eff), effect.mean, effect.sd) #what are the effect sizes?
  # (x - min(x))/(max(x) - min(x))
  return(eff)
}
## fixed number of effects, normal dist. For equal effects, can just do effect.sd = 0.
effect.dist.func.fn <- function(n){
  eff <- numeric(n)
  eff[sample(n, n.eff, replace = F)] <- rnorm(n.eff, effect.mean, effect.sd)
  return(eff)
}
## set the dist to use:
if(effect.dist == "zero.inflated.normal"){
  effect.dist.func <- effect.dist.func.znorm
}
if(effect.dist == "fixed.n.normal"){
  effect.dist.func <- effect.dist.func.fn
}
if(effect.dist == "bayesB"){
  effect.dist.func <- function(n) rbayesB(n, pi, d.f, scale)
}
#====================make data========================================
# read in data
x <- process_ms(x, chrl)
meta <- x$meta
x <- x$x

# assign effects
meta$effect <- rbayesB(nrow(meta), pi, d.f, scale)
cat("Mean effect size:", mean(meta$effect), "\n")

# assign phenotypes
phenos <- get.pheno.vals(x, meta$effect, .75)
saveRDS(list(x = x, meta = meta, phenos = phenos), outname)