args <- commandArgs(TRUE)
library(ggplot2)

# x <- as.character(args[1])
# outname <- as.character(args[2])
# run_n <- as.character(args[3])

# file and storage information:
runID <- "r01" # run ID and directory name where intermediate results will be stored. Will be created if needed.
x <- "data/theta4k_20000_1_rho40k.txt" # name of input ms format data
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

# trials for test
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


#========================run first prediction========================
#read in data
x <- process_ms(x, chrl)
meta <- x$meta
x <- x$x



#assign effects
# meta$effect <- effect.dist.func(nrow(meta))
# cat("Mean effect size:", mean(ch1.meta$effect), "\n")


ch1 <- x[meta$group == "chr1",]
ch1.meta <- meta[meta$group == "chr1",]


ch1.meta$effect <- effect.dist.func(nrow(ch1.meta))
cat("Mean effect size:", mean(ch1.meta$effect), "\n")

phenos <- get.pheno.vals(ch1, ch1.meta$effect, h)

# run prediction
pred_vals <- pred(x = as.matrix(ch1),
                  meta = ch1.meta,
                  phenotypes = phenos$p,
                  chr.length = chrl,
                  prediction.program = "ranger",
                  prediction.model = "RJ",
                  make.ig = make.ig,
                  sub.ig = F,
                  maf.filt = 0.05, 
                  runID = runID,
                  save.meta = save.meta, par = 10)

# refine once
imps <- pred_vals$output.model$model$variable.importance
b.imps <- quantile(abs(imps), .99)
keep <- which(abs(imps) >= b.imps)
x.keep <- pred_vals$x[keep,]
n.phenos <- phenos


pv_r <- pred(x = as.matrix(x.keep),
             meta = pred_vals$meta[keep,],
             phenotypes = n.phenos$p,
             chr.length = chrl,
             prediction.program = "ranger",
             prediction.model = "RJ",
             make.ig = make.ig,
             sub.ig = F,
             maf.filt = 0.05, 
             runID = runID,
             save.meta = save.meta, par = 10)


gens <- 100
vdat <- data.frame(gen = 1:gens, real_var = numeric(gens), pred_var1 = numeric(gens),
                   pred_var2 = numeric(gens), cor1 = numeric(gens), cor2 = numeric(gens))
test.x <- ch1
comp.dat <- array(dim = c(gens, ncol(ch1)/2, 3))
recom <- T

for(i in 1:gens){
  print(i)
  # do recomb
  if(i != 1){
    if(recom){
      test.x <- rand.mating(test.x, ncol(ch1)/2, meta = ch1.meta, 
                            rec.dist = rec.dist, chr.length = chrl, do.sexes = T)
    }
    else{
      test.x <- ch1[,sample(1:ncol(ch1), size = ncol(ch1), replace = T)]
    }
  }
  else{
    test.x <- ch1
  }
  
  # get real values
  t.phenos <- get.pheno.vals(test.x, ch1.meta$effect, h)
  vdat[i, "real_var"] <- var(t.phenos$a)
  
  # predict off of the first model
  t.dat <- convert_2_to_1_column(test.x[pred_vals$kept.snps,])
  colnames(t.dat) <- pred_vals$output.model$model$forest$independent.variable.names
  m1.phenos <- predict(pred_vals$output.model$model, data = t.dat)$predictions
  vdat[i, "pred_var1"] <- var(m1.phenos)
  vdat[i, "cor1"] <- cor(m1.phenos, t.phenos$a)
  
  # predict off of the second model
  t.dat <- convert_2_to_1_column(test.x[pred_vals$kept.snps,][keep,][pv_r$kept.snps,])
  colnames(t.dat) <- pv_r$output.model$model$forest$independent.variable.names
  m2.phenos <- predict(pv_r$output.model$model, data = t.dat)$predictions
  vdat[i, "pred_var2"] <- var(m2.phenos)
  vdat[i, "cor2"] <- cor(m2.phenos, t.phenos$a)
  
  # save bvs
  comp.dat[i,,] <- matrix(cbind(real = t.phenos$a, p1 = m1.phenos, p2 = m2.phenos))
}

vs <- vdat[,c("gen", "real_var", "pred_var1", "pred_var2")]
vs <- reshape2::melt(vs, id.vars = "gen")
colnames(vs)[2:3] <- c("model", "variance")
#ggplot(vs[vs$model != "real_var",], aes(x = gen, y = variance, color = model)) + geom_point()

cors <- vdat[,c("gen", "cor1", "cor2")]
cors <- reshape2::melt(cors, id.vars = "gen")
colnames(cors)[2:3] <- c("model", "correlation")
#ggplot(cors, aes(x = gen, y = correlation, color = model)) + geom_point()

save.image("rf_data_theta4k_20000_1_rho40k.RDS")