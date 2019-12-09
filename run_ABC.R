args <- commandArgs(TRUE)
x <- as.character(args[1])
#x <- "C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt"
outname <- as.character(args[2])
ABC_iter <- as.character(args[3])

# file and storage information:
save.meta <- F # should we save the metadata file when doing prediction/GWAS

# genome information
chrl <- 10000000 # chromosome length CURRENTLY MUST BE THE SAME ACROSS ALL CHRs

# effect size information
effect.dist <- "bayesB" # which effect distribution should we use
prob.effect <- 0.0001 # probability that any single SNP has an effect, for effect dists without a fixed number of effect loci.
h <- 0.75 # h^2, or heritability for the trait prior to any selection.
pi <- 1 - prob.effect # probability that a site has zero effect for "bayesB"
d.f <- 5 # degrees of freedom for the "bayesB" t-distribution
scale <- 1 # scale factor for the "bayesB" t-distribution

# chain information
chain_length <- 10000
burnin <- 1000
thin <- 100

# ABC params
## priors
pi_func <- function(x) rbeta(x, 200, 1)
df_func <- function(x) runif(x, 1, 100)
scale_func <- function(x) runif(x, 1, 1)

#========================set simulation parameters====================
# growth model information
r <- 2 # growth rate of population
K <- 500 # population carrying capacity

# survival model information
survival.dist <- "historic.variance"
fixed.var <- 45 # if a fixed survival variance (the "fixed.variance" model) is used, what should the survival varaince be?
sigma <- 1 # if the "historic.variance" option is used, the sigma (selection strength) variable of the gaussian kernal will be equal to sigma*historic_genetic_variance.

# selection optimum model
sopt.model <- "fixed" # which model? "fixed": fixed slide. "starting.variance" slide by a proportion of starting variance.
sopt.slide.iv <- 0.075 # by what proportion of initial genetic variance should the slection optimum slide each gen (for the "starting.variance" model) 
sopt.slide.fixed <- .6 # by how much should the survival option slide each gen (for the "fixed" model)
var.theta <- 0.1 # how much environmental stochasticity is there?

# simulation parameters
max_gens <- 100 # maximum number of gens to do selection
n_runs <- 10 # number of times to run the selection operation

#========================source scripts===============================
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.4", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(methods); library(data.table); library(doParallel); library(dplyr)

source("~/coalescence/prediction_accuracy/genomic_prediction_accuracy/growth_sim.R")


#========================models and distributions for simulations=======================================

# population growth
l_g_func <- function(x){
  return((K*x*exp(r))/(K + x*(exp(r) - 1)))
}

#surivival probability follows a gaussian distribution around the optimal phenotype.
## surivial variance based on historical genomic variance
s_gauss_scaled_func <- function(x, opt_pheno, hist.var = sqrt(h.pv), ...){
  x <- exp(-(x-opt_pheno)^2/(2*(hist.var*sigma)^2))
  return(x)
}
## fixed survival variance
s_gauss_fixed_var_scaled_func <- function(x, opt_pheno, ...){
  x <- exp(-(x-opt_pheno)^2/(2*fixed.var^2))
  return(x)
  #(x-min(x))/(max(x)-min(x)) #scaled.
}
## set appropriate model
if(survival.dist == "historic.variance"){
  survival.dist.func <- s_gauss_scaled_func
}
if(survival.dist == "fixed.variance"){
  survival.dist.func <- s_gauss_fixed_var_scaled_func
}

# selection shift
## increases the selection optimum by a percentage of starting variance
sopt_ivar_func <- function(x, iv,  slide = sopt.slide.iv, ...){ #increase is a percentage of starting variance
  if(iv == 0){stop("With a genetic variance of zero, the mean phenotype will not change each gen! Consider using a fixed phenotype slide.\n")}
  x <- x + iv*slide
}
## increases the selection optimum by a fixed amount
sopt_const_func <- function(x, slide = sopt.slide.fixed, ...){ #fixed increase, probably more realistic...
  x <- x + slide
}
## set the correct model
if(sopt.model == "fixed"){
  sopt.func <- sopt_const_func
}
if(sopt.model == "starting.variance"){
  sopt.func <- sopt_ivar_func
}

#recombination distribution, Poission with one expected recombination event per chr.
rec.dist <- function(x){
  return(rpois(x, lambda = 1))
}


#========================read in genomic data, assign effects, run ABC======

# 
# 
# # read in data
# x <- process_ms(x, chrl)
# meta <- x$meta
# x <- x$x
# 
# # assign effects
# meta$effect <- rbayesB(nrow(meta), pi, d.f, scale)
# cat("Mean effect size:", mean(meta$effect), "\n")
# 
# # assign phenotypes
# phenos <- get.pheno.vals(x, meta$effect, .75)
# saveRDS(list(x = x, meta = meta, phenos = phenos), "ABC/ABC_test_input_data.RDS")

d <- readRDS(x)
x <- d$x
meta <- d$meta
phenos <- d$phenos
rm(d); gc(); gc();

# run
ABC_res <- ABC_on_hyperparameters(x, phenos$p, iters = 1, pi_func = pi_func, 
                                  df_func = df_func, scale_func = scale_func, ABC_scheme = "C",
                                  h = h, chain_length = chain_length, burnin = burnin, thin = thin, par = F)

# save
ABC_res$dists$iter <- ABC_iter
write.table(ABC_res$dists, paste0("dists_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ABC_res$effects, paste0("effects_", ABC_iter, "_", outname), sep = "\t", quote = F, col.names = F, row.names = F)

#====================run simulations====================
# prep inputs
## real
sim.real.x <- list(e.eff = data.frame(site = 1:nrow(meta), effect = meta$effect), h = h,
                   phenotypes = phenos, meta = meta, x = x)
## pseudo
pseudo.meta <- meta
pseudo.meta$effect <- ABC_res$effects[[1]]
pseudo.phenos <- get.pheno.vals(x, pseudo.meta$effect, h)
sim.psuedo.x <- list(e.eff = data.frame(site = 1:nrow(meta), effect = pseudo.meta$effect), h = h,
                     phenotypes = pseudo.phenos, meta = pseudo.meta, x = x)

# prep outputs
out.gs.real <- vector("list", length = n_runs)
out.gs.pseudo <- out.gs.real

# run
for(i in 1:n_runs){
  # real
  out.gs.real[[i]] <- gs(x = sim.real.x,
                         gens = max_gens,
                         pred.method = "real",
                         growth.function =  l_g_func,
                         survival.function = survival.dist.func,
                         selection.shift.function = sopt.func,
                         rec.dist = rec.dist,
                         var.theta = var.theta,
                         plot_during_progress = F,
                         chr.length = chrl,
                         print.all.freqs = F,
                         fgen.pheno = T)$run_vars
  out.gs.real[[i]] <- cbind.data.frame(run = i, as.data.frame(out.gs.real[[i]]), source = "real", iter = ABC_iter, stringsAsFactors = F)
  
  # pseudo
  out.gs.pseudo[[i]]  <- gs(x = sim.psuedo.x,
                          gens = max_gens,
                          pred.method = "real",
                          growth.function =  l_g_func,
                          survival.function = survival.dist.func,
                          selection.shift.function = sopt.func,
                          rec.dist = rec.dist,
                          var.theta = var.theta,
                          plot_during_progress = F,
                          chr.length = chrl,
                          print.all.freqs = F,
                          fgen.pheno = T)$run_vars
  out.gs.pseudo[[i]] <- cbind.data.frame(run = i, as.data.frame(out.gs.pseudo[[i]]), source = "pseudo", iter = ABC_iter, stringsAsFactors = F)
}

# bind
out.gs.real <- dplyr::bind_rows(out.gs.real)
out.gs.pseudo <- dplyr::bind_rows(out.gs.pseudo)

# save
write.table(out.gs.real, paste0("gs_real_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(out.gs.pseudo, paste0("gs_pseudo_", outname), sep = "\t", quote = F, col.names = F, row.names = F)
