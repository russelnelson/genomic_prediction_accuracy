# This script reads in a set of simulated genomic data, randomly assigns phenotypic effects, runs
# genomic prediction on those effects, and then simulates selection forward in time for both the predicted
# and actual effects.

# Note: JWAS is very slow for large numbers of markers. Try GBLUP method in JWAS?

#========================define parameters==========
args <- commandArgs(TRUE)
library(ggplot2)

# x <- as.character(args[1])
# outname <- as.character(args[2])
# run_n <- as.character(args[3])

runID <- "r01" # run ID and directory name where intermediate results will be stored. Will be created if needed.
x <- "theta4k_1000_10_rho40k.txt" # name of input ms format data
outname <- "trial_300_runs.RDS" # file name to save output dataset NOT CURRENTLY IMPLEMENTED
chrl <- 10000000 # chromosome length CURRENTLY THE SAME ACROSS ALL CHRs
max_gens <- 1000 # maximum number of gens to do selection
n_runs <- 3 # number of times to run the selection operation NOT CURRENTLY IMPLEMENTED
h <- 0.5 # h^2, or heritability for the trait prior to any selection.
prob.effect <- 0.01 # probability that any single SNP has an effect
effect.sd <- 0.5 # sd of effect sizes
effect.mean <- 0 # mean of effect sizes
r <- 2 # growth rate of population NOT CURRENTLY IMPLEMENTED, CHANGE IN MODEL BELOW
make.ig <- TRUE # should new input files for JWAS be created?
sub.ig <- FALSE # Numeric or FALSE. Should we randomly subset out some SNPs to run in JWAS?
maf <- FALSE # Numeric or FALSE. Should we remove SNPs with a low minor allele frequency? If so, how low?
qtl_only <- TRUE # Should we only put QTLs into JWAS? If TRUE, overrides sub.ig and maf.
method <- "BGLR" # Which method are we using to generate estimated effect sizes?
model <- "BayesC" # What model should we use?
chain_length <- 20000 # How long should the MCMC chain in JWAS be?
burnin <- 2000 # How many MCMC iterations should we discard at the start of the chain? Must be less than chain_length!
thin <- 300 # How should the MCMC iterations be thinned?
pass.resid <- 0.25 # Numeric or NULL. Should we pass residual variance to JWAS? If a number, by what maximum proportion should we randomly adjust the value before passing it?
pass.var <- 0.25 # Numeric or NULL. Should we pass genetic variance to JWAS? If a number, by what maximum proportion should we randomly adjust the value before passing it?
standardize <- FALSE # Should phenotypes be centered and scaled before being passed to JWAS?
adjust_phenotypes <- TRUE # When doing selection on the prediced effect sizes, do we need to re-adjust phenotypic variance to match the initial variance? Mostly useful since variances predicted by JWAS are less than that in the data provided to it, so the variance needs to be adjusted if you want to use the observed phenotypes for the first round of selection.
intercept_adjust <- T # For selection on the predicted effect sizes, should we re-adjust phenotypes given an intercept (the mean phenotypic value of the input data)? 
julia.path <- "/Users/Hemstrom/AppData/Local/Julia-0.7.0/bin/julia.exe" # What is the path to julia.exe?
# path_libjulia <- "C:/Users/Hemstrom/AppData/Local/Julia-0.7.0/bin/libjulia.dll"
# JWASr::jwasr_setup_win(path_libjulia)

#========================models and distributions for simulations=======================================
#loci effect size distribution
effect.dist.func <- function(n){
  x <- rnorm(n, effect.mean, effect.sd)
  # (x - min(x))/(max(x) - min(x))
}

#logistic growth
l_g_func <- function(x, K =575, r = 2){
  return((K*x*exp(r))/(K + x*(exp(r) - 1)))
}

#surivival probability follows a normal distribution around the optimal phenotype.
s_norm_scaled_func <- function(x, opt_pheno, hist_var = h.pv, ...){
  x <- dnorm(x, opt_pheno, sqrt(hist_var)*4) #normal dist, sd is based on ancestral var
  ((.7-0)*(x-min(x))/(max(x) - min(x))) #scaled between 0 and .7
  #(x-min(x))/(max(x)-min(x)) #scaled.
}
s_norm_fixed_var_scaled_func <- function(x, opt_pheno, fixed_var = 1, ...){
  x <- dnorm(x, opt_pheno, fixed_var) #normal dist, sd is provided. 
  ((.7-0)*(x-min(x))/(max(x) - min(x))) #scaled between 0 and .5
  #(x-min(x))/(max(x)-min(x)) #scaled.
}


#increases the selection optimum by a percentage each gen
sopt_sp_func <- function(x, iv, slide = .1, ...){ #increase is a percentage of starting variance
  if(iv == 0){stop("With a genetic variance of zero, the mean phenotype will not change each gen! Consider using a fixed phenotype slide.\n")}
  x <- x + iv*slide
}
sopt_const_func <- function(x, slide = 1, ...){ #fixed increase, probably more realistic...
  x <- x + slide
}

#recombination distribution, Poission with one expected recombination event per chr.
rec.dist <- function(x){
  return(rpois(x, lambda = 1))
}


#========================source scripts===============================
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.4", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(methods); library(data.table)

source("growth_sim.R")

#========================read in genomic data and assign effects======

#read in data
x <- process_ms(x, chrl)
meta <- x$meta
x <- x$x

#assign effects
meta$effect <- rbinom(nrow(meta), 1, prob.effect) #does each site have an effect?
meta$effect[which(meta$effect == 1)] <- effect.dist.func(sum(meta$effect)) #what are the effect sizes?
cat("Mean effect size:", mean(meta$effect), "\n")

#generate individual effects
ind.effects <- get.pheno.vals(x, meta$effect, h = h, standardize = standardize)

#=========================call the prediction function====================
pred_vals <- pred(x = x, 
                  effect.sizes = meta$effect,
                  ind.effects = ind.effects,
                  chr.length = chrl,
                  method = method,
                  chain_length = chain_length,
                  burnin = burnin,
                  thin = thin,
                  model = model,
                  make.ig = make.ig, 
                  sub.ig = sub.ig, 
                  maf.filt = maf, 
                  julia.path = julia.path, 
                  runID = runID, 
                  qtl_only = qtl_only,
                  pass.resid = pass.resid,
                  pass.var = pass.var,
                  standardize = standardize)


# pred_vals <- readRDS("full_data_pred_vals.RDS") 
# # diagnostic plots
# ggplot(pred_vals$e.eff, aes(V2)) + geom_histogram()
# comp <- cbind(pred_vals$e.eff, pred_vals$meta)
# colnames(comp)[1:2] <- c("tag", "predicted")
# mcomp <- reshape2::melt(comp, id.vars = c("group", "position", "tag"))
# ggplot(mcomp, aes(x = position, y = value, color = variable)) + geom_point() + facet_wrap(~group)
# 
# 
# t_comp <- comp[,c(3,4,2)]
# colnames(t_comp)[3] <- "effect"
# t_comp <- cbind(t_comp, type = "predicted")
# bcomp <- rbind(t_comp, cbind(meta, type = "observed"))
# sm_meta <- snpR::s_ave_multi(meta, "effect", 200, 150, FALSE, levs = "group")
# colnames(sm_meta)[3] <- "effect"
# sm_meta <- sm_meta[,-4]
# 
# bcomp <- rbind(t_comp, cbind(sm_meta, type = "observed"))
# ggplot(bcomp, aes(x = position, y = effect, color = type)) + geom_line() + facet_wrap(~group)
# 
# sm_comp <- snpR::s_ave_multi(t_comp, "effect", 200, 150, FALSE, levs = "group")
# colnames(sm_comp)[3] <- "effect"
# sm_comp <- sm_comp[,-4]
# bcomp2 <- rbind(cbind(sm_comp, type = "predicted"), cbind(sm_meta, type = "observed"))
# ggplot(bcomp2, aes(x = position, y = effect, color = type)) + geom_line() + facet_wrap(~group)







#=========================prepare and run simulations for both the real and predicted effects=============
# note, using phenotypes created in the pred function as the first gen phenotypes for both runs!
pred_vals$x <- as.data.table(pred_vals$x)
x <- as.data.table(x)

#simulation on predicted data. The first generation phenotypes aren't used here, since they are far from what would be predicted.
pr.gs <- gs(x = pred_vals$x, 
            effect.sizes = pred_vals$e.eff[,2], 
            h = pred_vals$h, 
            gens = max_gens,
            growth.function = l_g_func, 
            survival.function = s_norm_scaled_func, 
            selection.shift.function = sopt_const_func, 
            rec.dist = rec.dist,
            meta = pred_vals$meta, 
            plot_during_progress = F, 
            chr.length = chrl, 
            print.all.freqs = F,
            adjust_phenotypes = adjust_phenotypes,
            intercept_adjust = intercept_adjust,
            fgen.pheno = pred_vals$a.eff$p)

#simulation on 'real' data
re.gs <- gs(x = x, 
            effect.sizes = meta$effect, 
            h = h, 
            gens = max_gens, 
            growth.function =  l_g_func, 
            survival.function = s_norm_scaled_func, 
            selection.shift.function = sopt_const_func, 
            rec.dist = rec.dist, 
            meta = meta, 
            plot_during_progress = F, 
            chr.length = chrl, 
            print.all.freqs = F,
            fgen.pheno = pred_vals$a.eff$p)

#========================return the result===============================

#plot summary data
pdat <- as.data.frame(rbind(cbind(model = "obs", gen = 1:nrow(re.gs$summary), as.data.frame(re.gs$summary)),
              cbind(model = "pred", gen = 1:nrow(pr.gs$summary), as.data.frame(pr.gs$summary))))
ggplot(pdat, aes(x = gen, y = N, color = model)) + theme_bw() + geom_line()
ggplot(pdat, aes(x = gen, y = var_a)) + theme_bw() + geom_line() + facet_wrap(~model, ncol = 1, scale = "free_y")
ggplot(pdat, aes(x = gen, y = diff)) + theme_bw() + geom_line() + facet_wrap(~model, ncol = 1, scale = "free_y")
ggplot(pdat, aes(x = gen, y = mu_a)) + theme_bw() + geom_line() + facet_wrap(~model, ncol = 1, scale = "free_y")

#plot allele data
# afdat <-  melt(cbind(pr.gs$frequencies, snpID = 1:nrow(pr.gs$frequencies), model = "pred"), id.vars = c("group", "position", "effect", "snpID", "model"))
# afdat <- rbind(afdat,
#                melt(cbind(pr.gs$frequencies, snpID = 1:nrow(pr.gs$frequencies), model = "obs"), id.vars = c("group", "position", "effect", "snpID", "model")))
# 
# 
# 
# 
# pr.as <- pr.gs$frequencies
# final.gen.col <- ncol(pr.as)
# pr.as$diff <- (pr.as[,final.gen.col]-pr.as[,4])/pr.as[,4]
# pr.as$status <- ifelse(pr.as[,final.gen.col] == 1, "fixed",
#                        ifelse(pr.as[,final.gen.col] == 0, "lost",
#                               "segregating"))
# 
# 
# pr.as$s_effect <- snpR::s_ave_multi(pr.as, "effect", 50, NULL, FALSE, "group")$smoothed_effect
# ggplot(pr.as, aes(x = position, y = s_effect, color = status)) + theme_bw() + geom_point() + facet_wrap(~group) + scale_color_viridis_d()
# 
# m1 <- lm(diff ~ effect, pr.as)
# summary(m1)
# 
# 
# pr.as <- pr.gs$frequencies
# pr.as$diff <- pr.as[,4] - pr.as[,ncol(pr.as)]
# ggplot(pr.as, aes(x = effect, y = diff)) + theme_bw() + geom_point()
# 
# 
# colnames(pf1dat)[5:6] <- c("gen", "freq")
# 
# ggplot(pf1dat, aes(x = position, y = freq, color = gen)) + geom_line() + facet_wrap(~group)
# 
# ggplot(pf1dat, aes(x = gen, y = freq, group = snpID, color = effect)) + geom_line() + scale_color_viridis_c() + theme_bw()
# 












#=================debugging============
e.dist.func <- function(A1, hist.a.var, h){
  esd <- sqrt((hist.a.var/h)-hist.a.var) # re-arrangement of var(pheno) = var(G) + var(E) and h2 = var(G)/var(pheno)
  env.vals <- rnorm(length(A1), 0, esd)
  return(env.vals)
}

plot(meta$effect)
plot(meta$effect*(rowSums(x)/ncol(x)))

# plot the addative genetic effects
temp <- weighted.colSums(x,meta$effect)
plot(temp)
temp <- temp[seq(1, length(temp), by = 2)] + temp[seq(2, length(temp), by = 2)]
plot(temp)

# plot the addative genetic effects of the loci taken for prediction
test.x <- x[which(paste0(meta$group, meta$position) %in% paste0(pred_vals$meta$group , pred_vals$meta$position)),]
test.meta <- meta[which(paste0(meta$group, meta$position) %in% paste0(pred_vals$meta$group , pred_vals$meta$position)),]

temp1 <- weighted.colSums(test.x,test.meta$effect)
plot(temp1)
temp1 <- temp1[seq(1, length(temp1), by = 2)] + temp1[seq(2, length(temp1), by = 2)]
plot(temp1)

#plot the addative genetic effects based on predicted effect sizes
temp2 <- weighted.colSums(pred_vals$x, pred_vals$e.eff$V2)
temp2 <- temp2[seq(1, length(temp2), by = 2)] + temp2[seq(2, length(temp2), by = 2)]
plot(temp2)
plot(temp2 + e.dist.func(temp2, var(temp2), pred_vals$h)) #including environmental effects

#comparison of both
plot(temp, temp2)
plot(pred_vals$a.eff$p, temp2 + e.dist.func(temp2, var(temp2), pred_vals$h))

#lots of residual variance that is genetic in nature is left by the model!



# influence of different loci on pop level effect of genetic value
plot(abs(pred_vals$e.eff$V2*(rowSums(pred_vals$x)/ncol(pred_vals$x))))
plot(abs(meta$effect[meta$effect != 0]*rowSums(x[meta$effect != 0,])/ncol(x))) #many more high impact loci
plot(abs(pred_vals$e.eff$V2*(rowSums(pred_vals$x)/ncol(pred_vals$x))),
     abs(meta$effect[meta$effect != 0]*rowSums(x[meta$effect != 0,])/ncol(x)))

#1: alleleic variance drops MASSIVELY after gen one, causing the effective h to plummet.
# here's why: rearranging the chromosomes randomly consistantly causes a severe drop in var(a)
#predicted data:
temp <-get.pheno.vals(pred_vals$x, effect.sizes = pred_vals$e.eff$V2, h = pred_vals$h) # note, h doesn't matter for this funciton, since we're only looking at genetic effects!
r.x <- pred_vals$x[,sample(x = ncol(pred_vals$x), ncol(pred_vals$x), replace = F)]
temp.r <- get.pheno.vals(as.data.table(r.x), pred_vals$e.eff$V2, h = pred_vals$h)
var(temp.r$a)/var(temp$a)

#simulated data:
temp <- get.pheno.vals(x, meta$effect, h = 1)
r.x <- x[,sample(1:ncol(x), ncol(x), replace = F)]
temp.r <- get.pheno.vals(r.x, meta$effect, h = 1)
var(temp.r$a)/var(temp$a)

# my guess: The model assumes a smaller variance in effect sizes than the data actually has. In order to best
# match the observed data, it therefore positions effect sizes to maximize the variance of the data set, since
# individual effects are not as large. See density plots of effect sizes:
plot(density(pred_vals$e.eff$V2)) # predicted, much lower variance in effect sizes
plot(density(pred_vals$meta$effect)) # real, much higher variance in effect sizes.
plot(pred_vals$e.eff$V2 ~ pred_vals$meta$effect) # real data has a much larger spread. Most of the predicted values are much closer to zero!




re.gs <- gs(x = x, 
            effect.sizes = meta$effect, 
            h = h, 
            gens = max_gens, 
            growth.function =  l_g_func, 
            survival.function = s_norm_scaled_func, 
            selection.shift.function = sopt_const_func, 
            rec.dist = rec.dist, 
            meta = meta, 
            plot_during_progress = F, 
            chr.length = chrl, 
            print.all.freqs = F)




pred_vals <- pred(x = x, 
                  effect.sizes = meta$effect,
                  ind.effects = ind.effects,
                  chr.length = chrl,
                  method = "BGLR",
                  chain_length = 20000,
                  burnin = 2000,
                  model = "BL",
                  make.ig = make.ig, 
                  sub.ig = sub.ig, 
                  maf.filt = maf, 
                  julia.path = julia.path, 
                  runID = runID, 
                  qtl_only = qtl_only,
                  pass.resid = pass.resid,
                  pass.var = pass.var,
                  standardize = standardize)

plot(density(pred_vals$e.eff$V2))
plot(density(meta$effect[meta$effect != 0]))
