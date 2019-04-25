#========================define parameters==========
args <- commandArgs(TRUE)
library(ggplot2)

# x <- as.character(args[1])
# outname <- as.character(args[2])
# run_n <- as.character(args[3])

# file and storage information:
runID <- "r01" # run ID and directory name where intermediate results will be stored. Will be created if needed.
x <- "C:/Users/hemst/Documents/ss_file_access/theta4k_1000_10_rho40k.txt" # name of input ms format data
save.meta <- T # should we save the metadata file when doing prediction/GWAS

# genome information
chrl <- 10000000 # chromosome length CURRENTLY MUST BE THE SAME ACROSS ALL CHRs

# effect size information
effect.dist <- "fixed.n.normal" # which effect distribution should we use
prob.effect <- 0.01 # probability that any single SNP has an effect, for effect dists without a fixed number of effect loci.
effect.sd <- .5 # sd of effect sizes
effect.mean <- 0 # mean of effect sizes
n.eff <- 200 # number of SNPs with effects for the "fixed.n.normal" model.
h <- 0.5 # h^2, or heritability for the trait prior to any selection.

# effect size estimation
maf <- FALSE # Numeric or FALSE. Should we remove SNPs with a low minor allele frequency? If so, how low?
qtl_only <- FALSE # Should we only put QTLs into JWAS/BGLR? If TRUE, overrides sub.ig and maf.
chain_length <- 20000 # How long should the MCMC chain in JWAS/BGLR be?
burnin <- 2000 # How many MCMC iterations should we discard at the start of the chain for JWAS/BGLR? Must be less than chain_length!
thin <- 300 # How should the MCMC iterations be thinned? For BGLR only.
make.ig <- TRUE # should new input files for JWAS be created?

# trials for test
trials <- 1000
outfile <- "Tests/g_var_decline/g_var_decline_1chr_no_rec_test.RDS"
outplot <- "Tests/g_var_decline/g_var_decline_1chr_no_rec_test.pdf"

# simulation parameters
adjust_phenotypes <- F # When doing selection on the prediced effect sizes, do we need to re-adjust phenotypic variance to match the initial variance? Mostly useful since variances predicted by JWAS/BGLR are less than that in the data provided to it, so the variance needs to be adjusted if you want to use the observed phenotypes for the first round of selection.
intercept_adjust <- F # For selection on the predicted effect sizes, should we re-adjust phenotypes given an intercept (the mean phenotypic value of the input data)? 


#========================models and distributions for simulations=======================================
#loci effect size distribution
## zero inflated normal
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
meta$effect <- effect.dist.func(nrow(meta))
cat("Mean effect size:", mean(meta$effect), "\n")


ch1 <- x[meta$group == "chr1",]
ch1.meta <- meta[meta$group == "chr1",]

phenos <- get.pheno.vals(ch1, ch1.meta$effect, h)

# run prediction
pred_vals <- pred(x = as.matrix(ch1),
                  h = .5,
                  meta = ch1.meta,
                  phenotypes = phenos$p,
                  chr.length = chrl,
                  prediction.program = "BGLR",
                  effect.sizes = ch1.meta$effect,
                  prediction.model = "BayesB",
                  make.ig = make.ig,
                  sub.ig = F,
                  maf.filt = 0.05, 
                  runID = runID, 
                  chain_length = 5000, 
                  burnin = 1000, 
                  thin = 100,
                  standardize = T,
                  save.meta = save.meta)


# Effect comparison
## effect sizes
effect.comp <- data.frame(position = ch1.meta$position, 
                          actual = ch1.meta$effect, 
                          predicted = pred_vals$e.eff$V2,
                          importance = ch1.meta$effect*(rowSums(ch1)/ncol(ch1)))
effect.comp <- reshape2::melt(effect.comp, id.vars = "position")
colnames(effect.comp)[2:3] <- c("source", "effect")
ggplot(effect.comp, aes(x = position, y = log(effect^2), color = source)) + geom_point() + theme_bw()

## individual genetic variance
ae <- pred.BV.from.model(pred_vals$output.model$mod, ch1, pred.method = "model", model.source = "BGLR", h = pred_vals$h)
a.comp <- data.frame(individual = 1:ncol(ch1), ae = ae$a, ar = phenos$a)
a.comp <- reshape2::melt(a.comp, id.vars = "individual")
colnames(a.comp)[2:3] <- c("source", "a")
ggplot(a.comp, aes(x = a)) + geom_histogram() + facet_wrap(~source, ncol = 1) + theme_bw()
plot(ae$a, phenos$a)



# attempt to use PCA to adjust the invidual variance distribution to be closer to the real data
## do a PCA on the raw data
# PCA <- prcomp(t(ch1))
# res <- summary(PCA)
# plot(res$importance[2,])
# plot(PCA$rotation[,1])

mafs <- rowSums(ch1)/ncol(ch1)
mafs[mafs > 0.5] <- 1 - mafs[mafs > 0.5]
maf.ch1 <- ch1[mafs >= 0.05,]

scor <- cor(t(maf.ch1))
e <- RSpectra::eigs_sym(scor, ncol(maf.ch1))


## get the lumps and their effects
#new.loadings <- t(t(e$vectors)/rowSums(abs(e$vectors)))
new.loadings <- e$vectors
new.scores <- pred_vals$e.eff$V2 * new.loadings # the rotation is are the loadings/eigen vectors, and contain the imporatance of each loci to the 1000 lumps.
new.scores <- colSums(new.scores) # these are the "lumps" with their effects

## get "genotypes" for gcs at lumps
lump.genos <- matrix(0, length(new.scores), ncol(ch1))
for(i in 1:ncol(ch1)){
  lump.genos[,i] <- colSums(maf.ch1[,i] * new.loadings)
}
new.a <- colSums(new.scores * lump.genos)
old.a <- colSums(ch1.meta$effect * ch1)
pred.a <- colSums(pred_vals$e.eff$V2 * pred_vals$x)
plot(pred.a, new.a) # well correlated with the old a
plot(old.a, new.a) # well correlated with the true a


## function to lump genos and get a
lump_and_bv <- function(genos, scores, loadings){
  lump.genos <- data.table::as.data.table(matrix(0, length(scores), ncol(genos)))
  genos <- data.table::as.data.table(genos)
  for(i in 1:ncol(genos)){
    set(lump.genos, j = as.integer(i), value = colSums(unlist(genos[,..i]) * loadings))
  }
  a <- colSums(scores * lump.genos)
  BV <- a[seq(1, length(a), by = 2)] +
    a[seq(2, length(a), by = 2)]
  return(BV)
}




# over a few generations, compare the variance with recombination
trials <- 5

var_test <- matrix(NA, nrow = trials, ncol = 7)
colnames(var_test) <- c("gen", "var_BV", "var_eBV", "var_pcaBV",
                        "cor_BV_eBV", "cor_BV_pcaBV", "cor_eBV_pcaBV")
BVs <- array(NA, dim = c(trials, ncol(ch1)/2, 3))
var_test[,1] <- 1:trials
recom <- ch1

for(i in 1:trials){
  print(i)
  # true bv
  BVs[i,,1] <- get.pheno.vals(recom, ch1.meta$effect, h)$a
  var_test[i,2] <- var(BVs[i,,1])
  
  # keep the correct values
  maf.recom <- recom[pred_vals$kept.snps,]
  
  # estimated bv
  BVs[i,,2] <- pred.BV.from.model(pred_vals$output.model$mod, maf.recom, "model", "BGLR", pred_vals$h, h.av = "fgen")$a
  var_test[i,3] <- var(BVs[i,,2])
  var_test[i,5] <- cor(BVs[i,,1], BVs[i,,2])
  
  # pca bv
  
  BVs[i,,3] <- lump_and_bv(maf.recom, new.scores, new.loadings)
  var_test[i, 4] <- var(BVs[i,,3])
  var_test[i, 6] <- cor(BVs[i,,1], BVs[i,,3])
  var_test[i, 7] <- cor(BVs[i,,2], BVs[i,,3])
  
  if(i != nrow(var_test)){
    recom <- rand.mating(recom, ncol(recom)/2, ch1.meta, rec.dist = rec.dist,
                         chrl, T, facet = "group")
  }
  
  gc();gc();gc();gc();gc()
}

library(reshape2)
m_var_test <- melt(as.data.frame(var_test), id.vars = "gen")
colnames(m_var_test) <- c("gen", "model", "variance")
result_plot <- ggplot(m_var_test, aes(x = gen, y = variance, color = model)) + geom_line() + theme_bw() + scale_color_viridis_d(option = "A")
result_plot

outfile <- "./Tests/g_var_decline/g_var_decline_PCA_test.RDS"
outplot <- "./Tests/g_var_decline/g_var_decline_PCA_test.pdf"
saveRDS(list(var_tes = var_test, BVs = BVs), outfile)
pdf(outplot)
result_plot
# Close the pdf file
dev.off() 
dev.off()

# actually worse than the vanilla approach!