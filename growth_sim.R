#=======internal functions=======
#function to generate random environmental effects for a given set of BVs, addative genetic variance (either historic or current from BVs are typical), and heritability.
# note that standardization isn't perfect, but will result in data with a mean very close to 0 and var close to 1. The randomness of assigning random environmental effects everywhere will make it imperfect
# downstream standardization of the phenotypes will fix this if using estimated effect sizes!
e.dist.func <- function(A1, hist.a.var, h, standardize = F){
  esd <- sqrt((hist.a.var/h)-hist.a.var) # re-arrangement of var(pheno) = var(G) + var(E) and h2 = var(G)/var(pheno)
  env.vals <- rnorm(length(A1), 0, esd)
  
  #if it standardization is requested, do so
  if(standardize){
    env.vals <- env.vals/sqrt(var(env.vals)/(1-h)) # set variance to 1 - h
    env.vals <- env.vals - mean(env.vals) # set mean to 0. Should be close, but not perfect because of the random draws.
  }
  return(env.vals)
}

#get phenotypic values given genotypes, effect sizes, and heritabilities. If hist.a.var is true, uses the amount of genomic variability this gen and h to figure out how big of an env effect to add. Otherwise uses the provided value (probably that in the first generation).
get.pheno.vals <- function(x, effect.sizes, h, hist.a.var = "fgen", standardize = FALSE){
  #get effect of each individual:
  a <- weighted.colSums(as.matrix(x), effect.sizes) # faster than t(x)%*%effect.sizes!
  
  a.ind <- a[seq(1, length(a), by = 2)] + a[seq(2, length(a), by = 2)] #add across both gene copies.
  
  #standardize the genetic variance if requested.
  if(standardize){
    a.ind <- a.ind/sqrt(var(a.ind)/h) # set the variance to h.
    a.ind <- a.ind - mean(a.ind) # set the mean to 0
  }
  
  #add environmental variance
  if(hist.a.var == "fgen"){
    pheno <- a.ind + e.dist.func(a.ind, var(a.ind), h, standardize)
    return(list(p = pheno, a = a.ind))
  }
  else{
    pheno <- a.ind + e.dist.func(a.ind, hist.a.var, h, standardize)
    return(list(p = pheno, a = a.ind))
  }
}

#converts 2 column to 1 column genotypes and transposes
convert_2_to_1_column <- function(x){
  if(!is.matrix(x)){x <- as.matrix(x)}
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  ind.genos <- matrix(ind.genos, nrow = ncol(x)/2, byrow = T) # rematrix and transpose!
  return(ind.genos)
}

#function to do column sums faster
src <- '
  Rcpp::NumericMatrix dataR(data);
  Rcpp::NumericVector weightsR(weights);
  int ncol = dataR.ncol();
  Rcpp::NumericVector sumR(ncol);
  for (int col = 0; col<ncol; col++){
  sumR[col] = Rcpp::sum(dataR( _, col)*weightsR);
  }
  return Rcpp::wrap(sumR);'

weighted.colSums <- inline::cxxfunction(
  signature(data="numeric", weights="numeric"), src, plugin="Rcpp")

# Function to calculate the estimated time untill a population begins to crash (growth rate less than one) based on Burger and Lynch 1995.
#    g_var: addative genetic variance
#    e_var: environmental variance
#    omega: width of the fitness function, usually given as omega^2
#    k: rate of environmental change in phenotypic standard deviations
#    B: mean number of offspring per individual
#    Ne: effective population size
#    theta_var: environmental stochasticity
B_L_t1_func <- function(g_var, e_var, omega, k, B, Ne, theta_var){
  # calc Vs
  Vs <- (omega^2) + e_var
  
  # calc Vlam
  # simplified: Vlam = (Vs*(1+2*Ne))/2*Ne + (((1+2*Vs)*(g_var+theta_var))/2*Vs)
  V_gt <- (Vs/(2*Ne)) + (g_var*theta_var)/(2*Vs)
  Vlam <- Vs + g_var + V_gt + theta_var
  
  #calc kc
  Bo <- B*omega/sqrt(Vlam)
  if(Bo < 1){
    return(list(t1 = NA, kc = NA, Vs = Vs, Vlam = Vlam, Bo = Bo))
  }
  kc <- (g_var/(g_var + Vs))*sqrt(2*Vs*log(Bo))
  
  if(k<kc){
    t1 <- Inf
  }
  else{
    t1 <- -((g_var + Vs)/g_var)*log(1-(kc/k))
  }
  
  #calc t1
  return(list(t1 = t1, kc = kc, Vs = Vs, Vlam = Vlam, Bo = Bo))
}


# function to extract BV predictions from a model
# pred.model: The GWAS/GP/ECT provided model
# g: the genotype matrix, where columns are gene copies
# pred.method: Should the BVs be predicted directly off the model ("model") or off of loci effects ("effects")?
# model.source: Which program made the model? Currently supports BGLR, JWAS, or RJ (ranger).
# h: heritability estimate
# h.av: historic genetic varaince, for prediction from effect sizes.
# effect.sizes: marker effect sizes, for prediction from effect sizes.
pred.BV.from.model <- function(pred.model, g, pred.method = NULL, model.source = NULL, h = NULL, h.av = "fgen", effect.sizes = NULL){
  if(pred.method == "effects"){
    pheno <- get.pheno.vals(g, effect.sizes, h, hist.a.var = h.av)
    a <- pheno$a #addative genetic values
    pheno <- pheno$p #phenotypic values
    return(list(a = a, p = pheno))
  }
  
  else{
    if(!is.matrix(g)){g <- as.matrix(g)}
    if(model.source == "BGLR"){
      g <- convert_2_to_1_column(g)
      a <- as.vector(g%*%pred.model$ETA[[1]]$b)
    }
    else if(model.source == "RJ"){
      g <- convert_2_to_1_column(g)
      colnames(g) <- pred.model$forest$independent.variable.names
      a <- predict(pred.model, as.data.frame(g))
    }
    else if(model.souce == "JWAS"){
      #in progress
    }
    else{
      #in progress
    }
    
    if(h.av == "fgen"){
      h.av <- var(a)
    }
    pheno <- a + e.dist.func(a, h.av, h)
    
    return(list(a = a, p = pheno))
  }
}


#=======distribution functions=========
#' Get random draws from the distribution used for bayesB regressions.
#'
#' Generate any number of random values drawn from the distribution used for
#' bayesB genomic regressions.
#'
#' Under a bayesB model, the effect of each site is drawn from a distribution
#' where var(g) == 0 with prob \emph{pi}, and is otherwise drawn from a scaled
#' t distribution with degrees of freedom \emph{d.f} and scale
#' \emph{scale}.
#'
#' @param n numeric. Number of draws to make from the distribution
#' @param pi numeric. Probability that any one site has zero effect
#' @param d.f numeric. Degrees of freedom for the scaled-inverse chi-squared
#'   distribution
#' @param scale numeric. Scale/shape parameter for the scaled-inverse
#'   chi-squared distribution.
rbayesB <- function(n, pi, d.f, scale){
  effects <- rbinom(n, 1, 1 - pi) # these are non-zero
  #effects[effects != 0] <- LaplacesDemon::rinvchisq(sum(effects), d.f, scale) # inverse chi distribution alternative
  effects[effects != 0] <- scale * rt(sum(effects), d.f)
  return(effects)
}

#=======function to do a single generation of random mating===========
rand.mating <- function(x, N.next, meta, rec.dist, chr.length, do.sexes = TRUE, facet = "group"){
  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }
  #=========get parents and assign gcs for the next gen====
  #make a new x with individuals in next gen
  ##find parents
  if(do.sexes){ # if there are two sexes
    sex <- rbinom(ncol(x)/2, 1, 0.5) #what are the parent sexes?
    if(sum(sex) == length(sex) | sum(sex) == 0){ # if every individual is the same sex, the population dies.
      return(NULL)
    }
    mates <- matrix(0, nrow = N.next, ncol = 2) # initialize, p1 and p2 are columns
    mates[,1] <- which(sex == 1)[sample(sum(sex), nrow(mates), T)] #get parents of sex a
    mates[,2] <- which(sex == 0)[sample(length(sex) - sum(sex), nrow(mates), T)] #get parents of sex b
  }
  else{ # if there is only one sex
    mates <- matrix(sample(ncol(x)/2, N.next*2, T), ncol = 2) #p1 and p2 are columns
    selfings <- which(mates[,1] == mates[,2]) #any selfing?
    while(length(selfings) > 0){ #correct selfing
      mates[selfings,] <- sample(ncol(x)/2, length(selfings)*2, T) #get new parents
      selfings <- which(mates[,1] == mates[,2]) #any selfing remaining?
    }
  }
  
  # table with the showing the distribution of the number of offspring for each adult:
  # table(c(table(mates), rep(0, (ncol(x)/2) - length(table(mates)))))
  x.next <- data.table::as.data.table(matrix(0, nrow(x), nrow(mates)*2)) #initialize x matrix for next gen
  
  #=========figure out which copy from each parent goes to offspring=====
  #randomly choose gene copies to push to individuals in the next gen. 
  # for each individual, do they get copy 1 or copy 2 from the parent?
  uf <- unique(meta[,facet])
  chr.source.p1.i <- matrix(rbinom(nrow(mates)*length(uf), 1, .5), ncol = length(uf), byrow = T) + 1
  chr.source.p2.i <- matrix(rbinom(nrow(mates)*length(uf), 1, .5), ncol = length(uf), byrow = T) + 1
  
  # add the correct copies
  ##which column does the data come from?
  chr.source.p1 <- ifelse(chr.source.p1.i == 1, mates[,1] * 2 - 1,
                          mates[,1] * 2)
  
  chr.source.p2 <- ifelse(chr.source.p2.i == 1, mates[,2] * 2 - 1,
                          mates[,2] * 2)
  
  #the other copies?
  chr.nsource.p1 <- ifelse(chr.source.p1.i == 1, mates[,1] * 2, 
                           mates[,1] * 2 - 1)
  chr.nsource.p2 <- ifelse(chr.source.p2.i == 1, mates[,2] * 2,
                           mates[,2] * 2 - 1)
  
  rm(chr.source.p1.i, chr.source.p2.i)
  
  #these now say which column in x to take the data from for each chr for each individual for bases that didn't recombine.
  
  #=========recombination and chromosome assignment======
  num.rec <- rec.dist(nrow(mates)*2*length(uf))
  rs <- sum(num.rec) #total number
  rec.pos <- sample(chr.length, rs, T) #get positions for each of these events. Assumes equal chr length, otherwise would need to put this into the per. chr loop. Could do later if necissary, but it'd be a bit slower.
  
  prog <- 0 #progress through num.rec tracker.
  pos.prog <- 0 #progress through recombination events tracker.
  
  #fill for each facet
  for(j in 1:length(unique(uf))){
    # cat(j,"\n")
    #overall approach: for each parent:
    # line up copy 1 and copy 2 chromosomes in a matrix, each column is a seperate chr. Copy one is the copy that is getting passed! Copy two is the one that isn't.
    # for each recombination event, on each chr, flip which chr we are taking from. For portions where we are taking from chr2, paste into chr 1, which will be the output.
    this.chr.pos <- meta$position[meta[,facet] == uf[j]]
    
    for(k in 1:2){
      
      trec <- c(prog + 1, prog + nrow(mates)) #which recombination events are we working with here?
      prog <- prog + nrow(mates)
      
      #get the number of recombination events per chr in this set
      tnr <- num.rec[trec[1]:trec[2]]
      
      #initialize matrix
      c1.mat <- data.table::as.data.table(matrix(0, length(this.chr.pos), ncol = nrow(mates)))
      c2.mat <- data.table::as.data.table(matrix(0, length(this.chr.pos), ncol = nrow(mates)))
      
      #paste in the values from x. c1.mat contains the "passed" chr, c2.mat contains the "unpassed" chr.
      if(k == 1){
        data.table::set(c1.mat, 1:nrow(c1.mat), as.integer(1:ncol(c1.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.source.p1[,j]])
        data.table::set(c2.mat, 1:nrow(c2.mat), as.integer(1:ncol(c2.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.nsource.p1[,j]])
      }
      else{
        data.table::set(c1.mat, 1:nrow(c1.mat), as.integer(1:ncol(c1.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.source.p2[,j]])
        data.table::set(c2.mat, 1:nrow(c2.mat), as.integer(1:ncol(c2.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.nsource.p2[,j]])
      }
      
      #only the recombining entries.
      wnz <- which(tnr != 0)
      if(length(wnz) == 0){
        #no recombination, mostly if pop size is VERY small...
        if(k == 1){
          data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(1, ncol(x.next), by = 2)), c1.mat)
        }
        else{
          data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(2, ncol(x.next), by = 2)), c1.mat)
        }
        next()
      }
      tnr_nz <- tnr[wnz]
      
      
      
      #get the positions of the recombination events
      trpos <- rec.pos[(pos.prog + 1):(pos.prog + sum(tnr_nz))]
      pos.prog <- pos.prog + sum(tnr_nz)
      
      
      # Now need to make and assign the actual output vector. I can't think of a good way to do this in an actually vectorized way, mostly because I have to look up positions to reference against for each chr.
      sort.prog.trpos <- 0 #how many positions have we searched through?
      for(m in 1:length(tnr_nz)){
        # cat("\t\t", m, "\n")
        sort.pos <- 
          c(sort(trpos[(sort.prog.trpos + 1):(sort.prog.trpos + tnr_nz[m])])) #sort the correct recomb positions and add the ending positions.
        sort.prog.trpos <- sort.prog.trpos + tnr_nz[m] #update.
        
        #now figure out which chr each position will be drawn from.
        chr.ident <- numeric(length(this.chr.pos))
        for(q in 1:length(sort.pos)){ 
          chr.ident <- chr.ident + (this.chr.pos <= sort.pos[q])
        }
        chr.ident <- chr.ident %% 2 #if the number of crossing over events was even, 0s mean copy 1 and 1s mean copy 2. Otherwise reversed.
        
        #assign.
        if(length(sort.pos) %% 2 == 0){
          data.table::set(c1.mat, which(chr.ident == 1), j = wnz[m], value = c2.mat[which(chr.ident == 1), wnz[m], with = FALSE])
          #assign entries where the chr flips in c1 mat to the respecitve entries in c2 mat.
        }
        else{
          data.table::set(c1.mat, which(chr.ident == 0), j = wnz[m], value = c2.mat[which(chr.ident == 0), wnz[m], with = FALSE])
        }
        
      }
      
      #================assign chromosomes to x.next.
      if(k == 1){
        data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(1, ncol(x.next), by = 2)), c1.mat)
      }
      else{
        data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(2, ncol(x.next), by = 2)), c1.mat)
      }
    }
  }
  
  return(x.next)
}

#=======function to do growth and selection=======
# note: init means are we simply intiating a population under selection. We'll need to keep all markers if so.
gs <- function(x, 
               gens, 
               growth.function, 
               survival.function, 
               selection.shift.function, 
               rec.dist,
               var.theta = 0,
               pred.method = "effects",
               plot_during_progress = FALSE, 
               facet = "group", chr.length = 10000000,
               fgen.pheno = FALSE,
               intercept_adjust = FALSE,
               print.all.freqs = FALSE,
               adjust_phenotypes = FALSE,
               do.sexes = TRUE,
               init = F){
  cat("Initializing...\n")
  #unpack x:
  if(pred.method == "effects"){ #unpack estimated effect sizes if provided.
    effect.sizes <- x$e.eff[,2]
  }
  if(fgen.pheno){ #unpack phenotypes if requested
    fgen.pheno <- x$phenotypes$p
  }
  h <- x$h
  meta <- x$meta
  if(pred.method != "real"){
    pred.mod <- x$output.model$mod
    pred.dat <- x$output.model$data
    model <- x$prediction.program
  }
  else{
    model <- "real"
    pred.method <- "effects" #since everything else works the same, just need to change inputs.
    effect.sizes <-  meta$effect
  }
  if(pred.method == "effects"){
    pred.mod <- NULL
    pred.dat <- NULL
  }
  if(pred.method == "model"){
    effect.sizes <- NULL
  }
  x <- x$x
  
  
  #=================checks========
  if(!pred.method %in% c("model", "effects")){
    stop("pred.method must be provided. Options:\n\tmodel: predict phenotypes directly from the model provided.\n\teffects: predict phenotypes from estimated effect sizes.\n")
  }
  
  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }
  
  if(pred.method == "effects"){
    if(nrow(x) != length(effect.sizes) | nrow(x) != nrow(meta)){
      stop("Provided x, effect sizes, and meta must all be of equal length!")
    }
  }
  else{
    if(nrow(x) != nrow(meta)){
      stop("Provided x and meta must be of equal length!")
    }
  }

  if(pred.method == "model"){
    if(!model %in% c("JWAS", "BGLR", "ranger")){
      stop("To predict from the model, a JWAS, BGLR, or ranger model must be provided.\n")
    }
  }
  else{
    if(model == "ranger"){
      stop("RF does not estimate effect sizes, so prediction must be done using the ranger model.\n")
    }
  }
  
  # before doing anything else, go ahead and remove any loci from those provided with no effect! Faster this way.
  # don't do this if initializing the population!
  if(pred.method == "effects" & !init){
    if(any(effect.sizes == 0)){
      n.eff <- which(effect.sizes == 0)
      x <- x[-n.eff,]
      meta <- meta[-n.eff,]
      effect.sizes <- effect.sizes[-n.eff]
    }
  }
  
  #=================get starting phenotypic values and BVs=========
  # get starting phenotypes and addative genetic values
  ## If first gen phenos aren't provided (should be uncommon)
  if(length(fgen.pheno) != ncol(x)/2){
    if(pred.method == "effects"){
      cat("Generating representative starting phenotypes from effect sizes.")
      pheno <- get.pheno.vals(x, effect.sizes, h)
      a <- pheno$a # BVs
      pheno <- pheno$p # phenotypic values
    }
    else{
      cat("Generating representative starting phenotypes from model.")
      a <- pred.BV.from.model(pred.mod, x, pred.method, model)
      pheno <- a +  e.dist.func(a, var(a), h) #add environmental effects
      #working here
    }
  }
  
  # otherwise use those, but still need to estimate BVs
  else{
    cat("Using provided phenotypic values.")
    pheno <- fgen.pheno #provded phenotypic values.
    
    a <- pred.BV.from.model(pred.model = pred.mod, g = x, pred.method = pred.method, 
                            model.source = model, h = h, h.av = "fgen", effect.sizes = effect.sizes)$a
    
  }
  
  #================set up BV variation adjustment to correct for drop in variance from GP methods============
  if(adjust_phenotypes){
    reorg_gcs <- rand.mating(x, ncol(x)/2, meta, rec.dist, chr.length, do.sexes, facet) # reorganize chrs once, since this causes one heck of a drop in var(a) in some GP results
    reorg_gcs <- rand.mating(reorg_gcs, ncol(x)/2, meta, rec.dist, chr.length, do.sexes, facet)
    re_p <- pred.BV.from.model(pred.model = pred.mod, g = reorg_gcs, pred.method = pred.method,
                                    model.source = model, h = h, h.av = "fgen", effect.sizes = effect.sizes) # re-predict BVs
    re_a <- re_p$a
    re_p <- re_p$p
    adj.a.var <- var(re_a) #variance next gen
    
    
    # now need to adjust a and pheno to fit the variance a gen later
    # multiply future phenotypes by the square root of these values, then adjust the mean back to the correct mean.
    
    
    
    # old version which adjusts back to the starting phenotypic var every generation.
    # reorg_gcs <- rand.mating(x, ncol(x)/2, meta, rec.dist, chr.length, do.sexes, facet) # reorganize chrs once, since this causes one heck of a drop in var(a) in some GP results
    # re_p <- pred.BV.from.model(pred.model = pred.mod, g = reorg_gcs, pred.method = pred.method, 
    #                                 model.source = model, h = h, h.av = "fgen", effect.sizes = effect.sizes) # re-predict BVs
    # re_a <- re_p$a
    # re_p <- re_p$p
    # ad.factor <- var(pheno)/(var(re_a)/h) # here's our adjustment factor
    # rm(re_a, re_p, reorg_gcs)
  }

  #if requested, get the amount to adjust phenotypes by in future gens.
  if(intercept_adjust){
    i.adj <- mean(pheno)
  }

  #================print out initial conditions, intiallize final steps, and run===========
  #starting optimal phenotype, which is the starting mean addative genetic value.
  #browser()
  opt <- mean(a) #optimum phenotype
  
  cat("\n\n===============done===============\n\nStarting parms:\n\tstarting optimum phenotype:", opt, 
      "\n\tmean phenotypic value:", mean(pheno), "\n\taddative genetic variance:", var(a), "\n\tphenotypic variance:", var(pheno), "\n\th:", h, "\n")
  
  #make output matrix and get initial conditions
  out <- matrix(NA, nrow = gens + 1, ncol = 8)
  colnames(out) <- c("N", "mu_pheno", "mu_a", "opt", "diff", "var_a", "stochastic_opt", "gen")
  N <- ncol(x)/2 #initial pop size
  h.av <- var(a) #get the historic addative genetic variance.
  h.pv <- var(pheno) #historic phenotypic variance.

  out[1,] <- c(N, mean(pheno), mean(a), opt, 0, h.av, opt, 1) #add this and the mean initial additive genetic variance
  if(plot_during_progress){
    library(ggplot2)
    pdat <- reshape2::melt(out)
    colnames(pdat) <- c("Generation", "var", "val")
    ranges <- data.frame(var = c("N", "mu_pheno", "mu_a", "opt", "diff"), 
                         ymin = c(0, out[1,2]*2, out[1,3]*2, out[1,4]*2, -10),
                         ymax = c(out[1,1]*1.05, 0, 0, 0, 10))
    pdat <- merge(pdat, ranges, by = "var")
    print(ggplot(pdat, aes(Generation, val)) + geom_point(na.rm = T) +
      facet_wrap(~var, ncol = 1, scales = "free_y", strip.position = "left") + 
      geom_blank(aes(y = ymin)) +
      geom_blank(aes(y = ymax)) +
      theme_bw() + xlim(c(0, max(pdat$Generation))) +
      theme(strip.placement = "outside", axis.title.y = element_blank(), strip.background = element_blank(),
            strip.text = element_text(size = 11)))
  }
  
  #initialize matrix to return allele frequencies if requested.
  if(print.all.freqs){
    a.fqs <- matrix(0, nrow(meta), gens + 1)
    a.fqs[,1] <- rowSums(x)/ncol(x)
  }
  
  #================loop through each additional gen, doing selection, survival, and fisher sampling of survivors====
  
  cat("\nBeginning run...\n\n================================\n\n")

  for(i in 2:(gens+1)){
    #=========survival====
    # get the optimum phenotype this gen
    t.opt <- rnorm(1, opt, var.theta)
    
    #survival:
    s <- rbinom(out[(i-1),1], 1, #survive or not? Number of draws is the pop size in prev gen, surival probabilities are determined by the phenotypic variance and optimal phenotype in this gen.
                survival.function(pheno, t.opt, hist.var = h.pv)) # calling the function in this way ensures that individuals with phenotypes at the optimum have a survival probability of whatever is set in the function.
    #if the population has died out, stop.
    if(sum(s) <= 1){
      if(print.all.freqs){
        a.fqs <- cbind(meta, a.fqs[,1:(i-1)], stringsAsFactors = F)
        out <- list(summary = out[1:(i-1),], frequencies = a.fqs)
      }
      else{
        out <- out[1:(i-1),]
      }
      return(list(run_vars = out, x = x, phenos = pheno, BVs = a))
    }
    
    #what is the pop size after growth?
    out[i,1] <- round(growth.function(sum(s)))
    
    #make a new x with the survivors
    x <- x[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors
    
    # # check phenotypic variance...
    # temp <- get.pheno.vals(x, effect.sizes, h, hist.a.var = h.av)
    # ptemp <- data.frame(val = c(a, temp$a), class = c(rep("T0", length(a)), rep("T1", length(temp$a))))
    # temp <- tem$p
    # if(intercept_adjust){
    #   temp <- temp + i.adj
    # }
    # # adjust variance
    # if(adjust_phenotypes != FALSE){
    #   s.p.mean <- mean(temp)
    #   temp <- temp*sqrt(ad.factor)
    #   temp <- temp - (mean(temp) - s.p.mean)
    # }
    # print(var(temp))
    
    #=============do random mating, adjust selection, get new phenotype scores, get ready for next gen====
    y <- rand.mating(x, out[i,1], meta, rec.dist, chr.length, do.sexes, facet)
    # check that the pop didn't die due to every individual being the same sex (rand.mating returns NULL in this case.)
    if(is.null(y)){
      return(list(run_vars = out, x = x, phenos = pheno, BVs = a))
    }
    else{
      x <- y
      rm(y)
    }
    
    #get phenotypic/genetic values
    pa <- pred.BV.from.model(pred.model = pred.mod, 
                             g = x, 
                             pred.method = pred.method, 
                             model.source = model, 
                             h = h, 
                             h.av = h.av, 
                             effect.sizes = effect.sizes)
    a <- pa$a
    pheno <- pa$p
    
    #if requested, adjust the phenotypic values.
    # adjust intercept
    if(intercept_adjust){
      pheno <- pheno + i.adj
    }
    # adjust variance
    if(adjust_phenotypes != FALSE){
      s.p.mean <- mean(pheno)
      pheno <- pheno*sqrt(ad.factor)
      pheno <- pheno - (mean(pheno) - s.p.mean)
    }
    
    #adjust selection optima
    opt <- selection.shift.function(opt, iv = sqrt(h.av))

    #save
    out[i,2] <- mean(pheno)
    out[i,3] <- mean(a)
    out[i,4] <- opt
    out[i,5] <- opt - mean(a)
    out[i,6] <- var(a)
    out[i,7] <- t.opt
    cat("gen:", i-1, 
        "\tf_opt:", round(out[i-1,4],3),
        "\ts_opt", round(out[i-1,7],3),
        "\tmean(pheno):", round(out[i,2],3),  
        "\tmean(a):", round(out[i,3],3),
        "\tvar(a):", round(var(a),3),
        "\tNs:", sum(s), 
        "\tN(t+1):", out[i,1],"\n")
    
    if(plot_during_progress){
      pdat <- reshape2::melt(out)
      colnames(pdat) <- c("Generation", "var", "val")
      ranges <- data.frame(var = c("N", "mu_pheno", "mu_a", "opt", "diff"), 
                           ymin = c(0, out[1,2]*2, out[1,3]*2, out[1,4]*2, -10),
                           ymax = c(out[1,1]*1.05, 0, 0, 0, 10))
      pdat <- merge(pdat, ranges, by = "var")
      print(ggplot(pdat, aes(Generation, val)) + geom_line(na.rm = T) + geom_point(na.rm = T) +
              facet_wrap(~var, ncol = 1, scales = "free_y", strip.position = "left") + 
              geom_blank(aes(y = ymin)) +
              geom_blank(aes(y = ymax)) +
              theme_bw() + xlim(c(0, max(pdat$Generation))) +
              theme(strip.placement = "outside", axis.title.y = element_blank(), strip.background = element_blank(),
                    strip.text = element_text(size = 11)))
    }
    
    
    #add allele frequencies if requested
    if(print.all.freqs){
      a.fqs[,i] <- rowSums(x)/ncol(x)
    }
    
    gc()
  }
  
  #prepare stuff to return
  browser()
  out[,"gen"] <- 1:nrow(out)

  if(print.all.freqs){
    a.fqs <- cbind(meta, a.fqs, stringsAsFactors = F)
    out <- list(summary = out, frequencies = a.fqs)
  }
  
  return(list(run_vars = out, x = x, phenos = pheno, BVs = a))
}


#=======function to process ms files=============
process_ms <- function(x, chr.length){
  infile <- x #infile
  lines <- readLines(x)
  lines <- lines[-which(lines == "")] #remove empty entries
  lines <- lines[-c(1,2)] #remove header info
  nss <- grep("segsites", lines) #get the number of segsites per chr
  chrls <- gsub("segsites: ", "", lines[nss]) #parse this to get the lengths
  chrls <- as.numeric(chrls)
  lines <- lines[-nss] #remove the segsites lines
  pos <- lines[grep("positions:", lines)] #find the positions
  lines <- lines[-grep("positions:", lines)] #remove the position
  div <- grep("//", lines) #find the seperators
  gc <- div[2] - div[1] - 1 #find the number of gene copies per chr
  if(is.na(gc)){gc <- length(lines) - 1} #if there's only one chr
  dat <- lines[-div] #get the data only
  dat <- strsplit(dat, "") #split the lines by individual snp calls
  x <- matrix(NA, nrow = sum(chrls), ncol = gc) #prepare output
  meta <- matrix(NA, nrow = sum(chrls), 2)
  
  #process this into workable data
  pchrls <- c(0, chrls)
  pchrls <- cumsum(pchrls)
  for(i in 1:length(chrls)){
    cat("\n\tChr ", i)
    tg <- dat[(gc*(i-1) + 1):(gc*i)] #get only this data
    tg <- unlist(tg) #unlist
    tg <- matrix(as.numeric(tg), ncol = chrls[i], nrow = gc, byrow = T) #put into a matrix
    tg <- t(tg) #transpose. rows are now snps, columns are gene copies
    tpos <- unlist(strsplit(pos[i], " ")) #grap and process the positions
    tpos <- tpos[-1]
    meta[(pchrls[i] + 1):pchrls[i + 1],] <- cbind(paste0(rep("chr", length = nrow(tg)), i), tpos)
    x[(pchrls[i] + 1):pchrls[i + 1],] <- tg #add data to output
  }
  
  meta <- as.data.frame(meta, stringsAsFactors = F)
  meta[,2] <- as.numeric(meta[,2])
  meta[,2] <- meta[,2] * chr.length
  
  colnames(meta) <- c("group", "position")
  colnames(x) <- paste0("gc_", 1:ncol(x))
  
  return(list(x = x, meta = meta))
}


#=======function to run gs in parallel and grab results======
pgs <- function(x, n_runs, effect.sizes, h, gens, growth.function, survival.function, 
                selection.shift.function, rec.dist,
                meta, facet = "group", chr.length = 10000000, par = "avail"){
  library(doParallel)
  
  if(par == "avail"){
    par <- parallel::detectCores(all.tests = TRUE)
  }
  
  cl <- snow::makeSOCKcluster(par, outfile="")
  doSNOW::registerDoSNOW(cl)
  
  ntasks <- n_runs
  progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
  opts <- list(progress=progress)
  
  output <- foreach::foreach(q = 1:ntasks, .inorder = FALSE,
                             .options.snow = opts, .export = "gs") %dopar% {
                               w_data <- gs(x, effect.sizes, h, gens, growth.function, survival.function, 
                                            selection.shift.function, rec.dist,
                                            meta, plot_during_progress = FALSE, facet, chr.length)
                             }
  
  parallel::stopCluster(cl)
  doSNOW::registerDoSNOW()
  
  res <- data.frame(gens_persisted = unlist(lapply(output, nrow)),
                    ending_a = unlist(lapply(output, function(x) x[nrow(x), 3])),
                    ending_optimum_a = unlist(lapply(output, function(x) x[nrow(x), 4])),
                    delta_a = unlist(lapply(output, function(x) abs(x[1,3] - x[nrow(x),3]))),
                    a_delta_rate = unlist(lapply(output, function(x) lm(mu_a ~ gen, cbind(gen = 1:nrow(x), as.data.frame(x)))$coefficients[2])),
                    end_var_a = unlist(lapply(output, function(x) x[nrow(x), 6])),
                    var_a_delta_rate = unlist(lapply(output, function(x) lm(var_a ~ gen, cbind(gen = 1:nrow(x), as.data.frame(x)))$coefficients[2])))
  
  names(output) <- paste0("run_", 1:length(output))
  output <- list(summary = res, results = output)
  
  return(output)
  
}



#=======function to initialize a population at a given selection optimum========
# wrapper for gs() with no shift in surival and assuming real effect sizes (no model). x needs to contain phenotypes and meta with effects
init_pop <- function(x,
                     init_gens, 
                     growth.function, 
                     survival.function, 
                     rec.dist,
                     var.theta = 0,
                     plot_during_progress = FALSE, 
                     facet = "group", chr.length = 10000000,
                     fgen.pheno = F,
                     print.all.freqs = FALSE,
                     do.sexes = TRUE){
  #=======set the parms for calling gs======
  s.shift.null <- function(x, ...){
    return(x)
  }
  
  if(!"effect" %in% colnames(x$meta)){
    stop("Effects must be provided in as a column in x$meta.\n")
  }
  if(!"h" %in% names(x)){
    stop("heritability must be provided in x$h.\n")
  }
  if(!is.numeric(x$h) | length(x$h) != 1){
    stop("heritability must be a single numeric value.\n")
  }
  #call gs correctly to initialize with no shift in the survival function and assuming real effects.
  out <- gs(x = x, 
            gens = init_gens, 
            growth.function = growth.function, 
            survival.function = survival.function, 
            selection.shift.function = s.shift.null, 
            rec.dist = rec.dist, 
            var.theta = var.theta, 
            pred.method = "real", 
            plot_during_progress = plot_during_progress, 
            facet = facet, 
            chr.length = chr.length, 
            fgen.pheno = fgen.pheno,
            intercept_adjust = F, 
            print.all.freqs = print.all.freqs, 
            adjust_phenotypes = F, 
            do.sexes = do.sexes,
            init = T)
  return(out)
}


#=======function to take in input data, use genomic prediction to estimate effect sizes based on phenotypes==============
#arugments:
#    x: Input data, matrix, df, or data.table. Converted internally to data.table. Columns are gene copies, rows are SNPs, formatted as 0 and 1.
#    effect.sizes: Vector of marker effect sizes. Will be optional eventually for prediction from phenotypes only.
#    ind.effects: Individual addative genetic values/BVs for individuals. Will eventually be optional for prediction from phenotypes only.
#    method: What method should we use to generate predictions? BGLR, JWAS, or ranger.
#    chain.length: How long should the MCMC chain in JWAS/BGLR be?
#    burnin: How many MCMC iterations should we discard at the start of the chain for JWAS/BGLR? Must be less than chain_length!
#    thin: How should the MCMC iterations be thinned? For BGLR only.
#    model: What model should we use? BGLR: FIXED, BL, BayesA-C, BRR. JWAS: G-BLUP, BayesA-C, RR-BLUP. ranger: RJ.
#    make.ig: should new input files for JWAS be created? If they don't exist already, set to TRUE.
#    sub.ig: 
#    pass.resid: NULL or numeric >= 0. A numeric value tells the function to pass the 
#                estimated residual variance in the model on to JWAS. A numeric value of 0 passes the exact variance,
#                a numeric value other than zero will fudge the variance number by up to the proportion given (1 fudges up to 100%).
#    pass.var: NULL or numeric >= 0. Like pass.resid, but for the true genetic variance.
#    standardize: Boolean. Should the addative genetic values be centered and scaled between -1 and 1 prior to entry into JWAS? Phenotypic values still won't be centered!
pred <- function(x, meta = NULL, effect.sizes = NULL, phenotypes = NULL,
                 prediction.program = "JWAS",
                 chain_length = 100000, 
                 burnin = 5000,
                 thin = 100,
                 prediction.model = NULL,
                 make.ig = TRUE, sub.ig = FALSE, maf.filt = 0.05, 
                 julia.path = "julia", runID = "r1", qtl_only = FALSE,
                 pass.resid = FALSE, pass.var = FALSE, 
                 ntree = 50000,
                 mtry = 1,
                 h = NULL,
                 standardize = FALSE,
                 save.meta = TRUE, par = NULL, pi = NULL, verbose = T){
  #============sanity checks================================
  # check that all of the required arguments are provided for the prediction.model we are running
  if(prediction.program %in% c("JWAS", "BGLR", "PLINK", "TASSEL", "ranger")){
    
    # JWAS checks
    if(prediction.program == "JWAS"){
      # chain length and burnin
      if(!is.numeric(c(chain_length, burnin))){
        stop("Chain_length and burnin must be numeric values.")
      }
      if(chain_length <= burnin){
        stop("Chain_length must be larger than burnin in order to estimate effect sizes.")
      }
      
      # path to julia
      if(!is.character(julia.path)){
        stop("Invalid path to julia executable.")
      }
      if(!file.exists(julia.path)){
        stop("Invalid path to julia executable.")
      }
      
      # JWAS prediction.model-for now, only RR-BLUP.
      if(!prediction.model %in% c("RR-BLUP", "BayesB")){
        stop("Invalid JWAS prediction.model.")
      }
      
      # check that there is prior info to pass if pass resid is ture
      if((pass.resid != FALSE | pass.var != FALSE) & !is.null(effect.sizes)){
        stop("Marker effect sizes must be defined in order to pass prior residual and variance info to JWAS.")
      }
      
    }
    
    # BGLR checks:
    else if(prediction.program == "BGLR"){
      # chain length and burnin
      if(!is.numeric(c(chain_length, burnin))){
        stop("Chain_length and burnin must be numeric values.")
      }
      
      # check for accepted prediction.model.
      if(!prediction.model %in% c("BRR", "FIXED", "BayesA", "BayesB", "BayesC", "BL")){
        stop(paste0("prediction.model ", prediction.model, " not recognized by BGLR. Options: BRR, FIXED, BayesA-C, BL.\n"))
      }
    }
    
    # random forest checks:
    else if(prediction.program == "ranger"){
      # prediction.model
      if(!(prediction.model %in% c("RJ"))){
        stop("For random forest, the prediction.model must be RJ (random jungle) for now.")
      }
      
      # ntree
      if(!is.numeric(ntree)){
        stop("ntree must be an integer!")
      }
      if(ntree != floor(ntree)){
        stop("ntree must be an integer!")
      }
      
      if(mtry > 1 | mtry < 0){
        stop("mtry must be between 1 and 0.\n")
      }
    }
    
  }
  else{
    stop("Invalid prediction.program provided. Options: JWAS, BGLR, PLINK, TASSEL, or ranger.")
  }
  
  
  # check that we are either provided with input marker effects or with input phenotypes
  if(all(c(is.null(phenotypes), is.null(effect.sizes)))){
    stop("Individual phenotypes or marker effect sizes must be provided!")
  }
  else if(all(c(is.null(phenotypes), is.null(effect.sizes)))){
    warning("Both phenotypes and marker effect sizes provided. Input phenotypes will be ignored!")
  }
  
  # check that qtl_only filtering isn't requested if effect sizes aren't provided!
  if(qtl_only & is.null(effect.sizes)){
    stop("qtl_only filtering can only be performed if marker effect sizes are provieded.")
  }
  
  if(!is.null(effect.sizes)){
    if(!is.numeric(h)){
      stop("h must be provided if effect.sizes are used.\n")
    }
  }
  
  if(is.null(meta[1]) & save.meta){
    save.meta <- F
    warning("Since no SNP metadata provided, no SNP metadata will be saved.\n")
  }
  
  cat("Preparing model inputs...\n")
  
  #============subfunctions=================================
  # do filters if requested
  filter_snps <- function(x, qtl_only, sub.ig, maf.filt, effect.sizes){
    rejects <- numeric(nrow(x)) # track rejected snps
    # qtl_only
    if(qtl_only){
      if(sub.ig != FALSE | maf.filt != FALSE){
        warning("No subsampling (sub.ig) or maf filtering (maf.filt) will occur when qtl_only = TRUE.\n")
      }
      s.markers <- which(effect.sizes != 0)
      x <- x[s.markers,]
      rejects[-s.markers] <- 1
    }
    
    # otherwise, other filters to check
    else{
      #filter low minor allele frequencies if requested (like many prediction studies will do!)
      if(maf.filt != FALSE){
        af <- matrixStats::rowSums2(x)/ncol(x)
        m.keep <- af >= maf.filt & af <= (1-maf.filt)
        x <- x[m.keep,]
        rejects[-which(m.keep)] <- 1
      }
      
      #subset markers
      if(sub.ig != FALSE){
        if(nrow(x) > sub.ig){
          s.markers <- sort(sample(nrow(x), sub.ig))
          x <- x[s.markers,]
          rejects[rejects != 1][-s.markers] <- 1
        }
        else{
          warning("Fewer markers than sub.ig. Running all markers.\n")
        }
      }
    }
    return(list(x = x, snp_ids = which(rejects == 0)))
  }
  
  
  #============prepare directories and phenotypes=========
  if(!is.null(phenotypes)){
    # if standardization is requested, set the mean phenoytpes to o and var(pheno) to 1
    if(standardize){
      phenotypes <- phenotypes/sd(phenotypes)
      phenotypes <- phenotypes - mean(phenotypes)
    }
    r.ind.effects <- list(p = phenotypes) # backup the ind effects.
  }
  
  # get phenotypes if effect sizes are provided.
  if(!is.null(effect.sizes)){
    phenotypes <- get.pheno.vals(x, effect.sizes, h = h, standardize = standardize)
    r.ind.effects <- phenotypes
    phenotypes <- phenotypes$p
  }
  
  # # create the directory to store results if it doesn't already exist and move over to it.
  if(!dir.exists(runID)){
    dir.create(runID)
  }
  owd <- getwd()
  setwd(runID)
  
  #============format data for prediction/GWAS==============
  #filter:
  x <- filter_snps(x, qtl_only, sub.ig, maf.filt, effect.sizes)
  kept.snps <- x$snp_ids
  if(!is.null(meta)){
    meta <- meta[kept.snps,]
  }
  x <- x$x
  
  
  if(prediction.program == "JWAS"){
    # make an individual effect file.
    ind.effects <- cbind(samp = as.character(1:500), phenotypes = phenotypes)
    write.table(ind.effects, "ie.txt", quote = F, col.names = T, row.names = F)
    
    # make an individual genotype file if it isn't already constructed.
    if(make.ig){
      # convert format
      ind.genos <- convert_2_to_1_column(x)
      ind.genos <- cbind(samp = 1:nrow(ind.genos), ind.genos) # add sample info
      colnames(ind.genos) <- c("samp", paste0("m", 1:(ncol(ind.genos)-1)))
      data.table::fwrite(ind.genos, "ig.txt", sep = " ", col.names = T)
    }
  }
  
  else if(prediction.program == "BGLR"){
    # convert
    ind.genos <- convert_2_to_1_column(x) # rows are individuals, columns are SNPs
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS
    
    # prepare ETA
    ETA <- list(list(X = ind.genos, model = prediction.model, saveEffects = T))
  }
  
  else if(prediction.program == "ranger"){
    t.x <- convert_2_to_1_column(x) # add sample info
    colnames(t.x) <- paste0("m", 1:ncol(t.x))
    
    t.eff <- data.frame(phenotype = phenotypes, stringsAsFactors = F)
    t.eff <- cbind(t.eff, t.x)
  }
  
  
  #=========run genomic prediction or GWAS and return results==========
  # for JWAS:
  if(prediction.program == "JWAS"){
    cat("Calling JWAS.\n")
    options(scipen = 999)
    julia.call <- paste0(julia.path, " ", owd, "/analysis.jl ", chain_length, " ", burnin)
    # add the residual and genetic variance if requested.
    if(!is.null(pass.resid)){
      rv <- var(r.ind.effects$p - r.ind.effects$a)
      rv <- rv + rv*runif(1, 0, pass.resid) #fudge according to factor provided
    }
    else{
      rv <- 1
    }
    if(!is.null(pass.var)){
      gv <- var(r.ind.effects$a)
      gv <- gv + gv*runif(1, 0, pass.var) #fudge according to factor provided
    }
    else{
      gv <- 1
    }
    julia.call <- paste0(julia.call, " ", rv, " ", gv, " ", prediction.model)
    if(!is.null(pi)){
      julia.call <- paste0(julia.call, " ", pi)
    }
    else{
      julia.call <- paste0(julia.call, " ", "false")
    }
    system(julia.call)
    
    #=========grab output and modify it to give the estimated effect size per locus=============
    e.eff <- read.table("est_effects.txt", header = F, sep = "\t")
    h <- read.table("h.txt")
    
    #save metadata for the selected markers if requested.
    if(save.meta){
      write.table(meta, "est_meta.txt", sep = "\t", quote = F, col.names = T, row.names = F)
    }
    
    setwd(owd)
    
    return(list(x = x, e.eff = e.eff, phenotypes = r.ind.effects, meta = meta, h = as.numeric(h), kept.snps = kept.snps))
  }
  
  # for BGLR:
  else if(prediction.program == "BGLR"){
    cat("Calling BGLR.\n")
    BGLR_mod <- BGLR::BGLR(y = phenotypes, ETA = ETA, nIter = chain_length, burnIn = burnin, thin = thin, verbose = verbose)
    
    # grab h2 estimate
    B <- BGLR::readBinMat('ETA_1_b.bin')
    h2 <- rep(NA,nrow(B))
    varU <- h2
    varE <- h2
    for(i in 1:length(h2)){
      u <- ind.genos%*%B[i,]	
      varU[i] <- var(u)
      varE[i] <- var(phenotypes-u)
      h2[i] <- varU[i]/(varU[i] + varE[i])
    }
    h2 <- mean(h2)
    
    # return the values
    e.eff <- data.frame(V1 = BGLR_mod$ETA[[1]]$colNames, V2 = BGLR_mod$ETA[[1]]$b, stringsAsFactors = F)
    write.table(e.eff, "est_effects.txt", sep = "\t", col.names = F, row.names = F, quote = F)
    write.table(h2, "h.txt", quote = F)
    if(save.meta){
      write.table(meta, "est_meta.txt", sep = "\t", quote = F, col.names = T, row.names = F)
    }
    
    setwd(owd)
    return(list(x = x, e.eff = e.eff, phenotypes = r.ind.effects, meta = meta, h = h2, prediction.program = "BGLR",
                prediction.model = prediction.model, output.model = list(mod = BGLR_mod, data = ETA), kept.snps = kept.snps))
  }
  
  # for randomForest
  else if(prediction.program == "ranger"){
    cat("Running randomforest (ranger rj implementaiton).\n")
    
    # figure out the mtry to use
    mtry <- nrow(x)*mtry
    
    # run the randomForest/jungle
    if(ncol(t.eff) - 1 >= 10000){
      rj <- ranger::ranger(dependent.variable.name = "phenotype", data = t.eff, mtry = mtry, 
                           num.trees = ntree, importance = "permutation", verbose = T, save.memory = T, num.threads = par)
    }
    else{
      rj <- ranger::ranger(dependent.variable.name = "phenotype", data = t.eff,
                           mtry = mtry, num.trees = ntree, importance = "permutation", verbose = T, num.threads = par)
    }
    
    return(list(x = x, phenotypes = r.ind.effects, meta = meta, prediction.program = "ranger",
                prediction.model = "RJ", output.model = list(model = rj), kept.snps = kept.snps))
  }
}



#=======function to do ABC on hyperparameters============
#' Conduct an ABC on a range of effect size distribution hyperparameters
#' 
#' Runs Approximate Bayesian Computation across a range of marker effect size
#' distribution hyperparameters using one of three different schemes in order to determine
#' the hyperparamters that generate a distribution most like the real genomic architecture of the trait.
#' 
#' ABC schemes: \itemize{
#'     \item{"A": }{Genomic Data -> prediction with prior hyperparameters -> compare predicted phenotypes to real phenotypes.}
#'     \item{"B": }{Genomic Data -> generate psuedo marker effects using distribution with prior hyperparameters -> prediction with defaults -> compare predicted phenotypes to real phenotypes.}
#'     \item{"C": }{Part 1: Genomic Data -> generate psuedo marker effects using distribution with prior hyperparameters -> "pseudo" predicted marker effects. 
#'                  Part 2: Genomic Data -> prediction with defaults -> direct estimated marker effects.
#'                  Part 3: Compare direct to "pseudo" estimated marker effects.}
#' }
#' 
#' @param x matrix. Input genotypes, SNPs as rows, columns as individuals. Genotypes formatted as 0,1,2 for the major homozygote, heterozygote, and minor homozygote, respectively.
#' @param phenotypes numeric vector. Observed phenotypes, one per individual.
#' @param iters numeric. Number of ABC permutations to run.
#' @param pi_func function, default function(x) rbeta(x, 25, 1). A distribution function for generating pi prior. Should take only one argument (n, the number of iters).
#' @param df_func function, default NULL. A distribution function for generating df prior. Should take only one argument (n, the number of iters).numeric vector of length 2, default NULL. Range (min, max) of degrees of freedom values to run.
#' @param scale_func function, default NULL. A distribution function for generating scale prior. Should take only one argument (n, the number of iters).numeric vector of length 2, default NULL. Range (min, max) of scale values to run.
#' @param h numeric, default NULL. Heritability to use. Will take a range in the future.
#' @param julia.path character, defualt "julia". File path to the julia executable, required for JWAS.
#' @param chain_length numeric, default 100000. Length of the MCMC chains used in each step of the ABC.
#' @param burnin numeric, default 5000. Number of MCMC chain steps discarded at the start of each MCMC.
#' @param thin numeric, default 100. Number of MCMC chain steps discarded between each sample used to form the posterior.
#' @param method character, default "bayesB". The marker effect size distribution/prediction method to use.
#' @param ABC_scheme character, default "A". The ABC_scheme to use. See details.
#' @param par numeric or FALSE, default FALSE. If numeric, the number of cores on which to run the ABC.
#' @param run_number numeric, default NULL. Controls how the itermediate output directories are named. If numeric, will be named for the number, otherwise, will be named for the iteration.
ABC_on_hyperparameters <- function(x, phenotypes, iters, pi_func = function(x) rbeta(x, 25, 1), 
                                   df_func = NULL,  scale_func = NULL, h = NULL,
                                   julia.path = "julia", chain_length = 100000, 
                                   burnin = 5000,
                                   thin = 100, method = "BayesB", ABC_scheme = "A", 
                                   par = F, run_number = NULL, est_h = F){

  
  # ks <- which(matrixStats::rowSums2(x)/ncol(x) >= 0.05)
  #============general subfunctions=========================
  euclid.dist <- function(p, o){
    dist <- sqrt((o - p)^2)
    return(dist)
  }
  euclid.distribution.dist <- function(p, o){
    dist <- ks.test(p, o)$statistic
    return(dist)
  }
  generate_pseudo_effects <- function(x, pi, df, scale, method, h = NULL){
    if(method == "BayesB"){
      pseudo_effects <- rbayesB(nrow(x), pi, df, scale)
      pseudo_phenos <- get.pheno.vals(x, pseudo_effects, h)$p
    }
    return(list(e = pseudo_effects, p = pseudo_phenos))
  }
  
  #============ABC_scheme functions for one rep=============
  scheme_A <- function(x, phenotypes, pi, method, t_iter){
    p <- pred(x, pi = pi, phenotypes = phenotypes, julia.path = julia.path, 
              burnin = burnin, thin = thin, chain_length = chain_length,
              prediction.program = "JWAS", prediction.model = method, runID = t_iter, verbose = F)
    
    dist <- euclid.dist(p$est.phenos, phenotypes)
    return(dist)
  }
  scheme_B <- function(x, phenotypes, pi, df, scale, method, t_iter){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h)
    
    p <- pred(x, phenotypes = pseudo$p,
              burnin = burnin, thin = thin, chain_length = chain_length,  
              prediction.program = "BGLR", prediction.model = method, runID = t_iter, verbose = F)
    
    p.phenos <- as.vector(convert_2_to_1_column(p$x)%*%p$output.model$mod$ETA[[1]]$b)
    
    dist <- euclid.distribution.dist(p.phenos, phenotypes)
    return(return(list(dist = dist, e = pseudo$e)))
  }
  scheme_C <- function(x, phenotypes, r.p.eff, pi, df, scale, method, t_iter){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h)
  
    cat("Beginning pseudo data", method, "run.\n")
    pseudo.pred <-pred(x, phenotypes = pseudo$p,
                       burnin = burnin, thin = thin, chain_length = chain_length,  
                       prediction.program = "BGLR", prediction.model = method,
                       runID = paste0(t_iter, "_pseudo"), verbose = F)
    
    dist <- euclid.distribution.dist(r.p.eff, pseudo.pred$output.model$mod$ETA[[1]]$b)
    return(list(dist = dist, e = pseudo$e))
  }
  scheme_D <- function(x, phenotypes, pi, df, scale, method){
    pseudo <- generate_pseudo_effects(x, pi, df, scale, method, h)
    dist <- euclid.distribution.dist(phenotypes, pseudo$p)
    return(list(dist = dist, e = pseudo$e))
  }
  
  loop_func <- function(x, phenotypes, pi, df, scale, method, scheme, t_iter, r.p.phenos = NULL, r.p.eff = NULL){
    if(scheme == "A"){
      dist <- scheme_A(x, phenotypes, pi, method, t_iter)
    }
    else if(scheme == "B"){
      dist <- scheme_B(x, phenotypes, pi, df, scale, method, t_iter)
    }
    else if(scheme == "C"){
      dist <- scheme_C(x, phenotypes, r.p.eff, pi, df, scale, method, t_iter)
    }
    else if(scheme == "D"){
      dist <- scheme_D(x, phenotypes, pi, df, scale, method)
    }
    return(dist)
  }
  
  #============ABC loop======================================
  # get the random values to run
  run_pis <- pi_func(iters)
  if(!is.null(df_func)){
    run_dfs <- df_func(iters)
  }
  else{
    run_dfs <- rep(NA, iters)
  }
  if(!is.null(scale_func)){
    run_scales <- scale_func(iters)
  }
  else{
    run_scales <- rep(NA, iters)
  }
  
  
  # initialize storage
  out <- cbind(pi = run_pis, df = run_dfs, scale = run_scales, dist = 0)
  
  # if doing a method where prediction needs to be run on the real data ONCE, or if h should be estimated, do that now:
  if(ABC_scheme == "C" | est_h == T){
    cat("Beginning real data", method, "run.\n")
    real.pred <- pred(x, phenotypes = phenotypes,
                      burnin = burnin, thin = thin, chain_length = chain_length,  
                      prediction.program = "BGLR", prediction.model = method, runID = "real_pred", verbose = F)
    if(ABC_scheme == "C"){
      #r.p.phenos <- as.vector(convert_2_to_1_column(real.pred$x)%*%real.pred$output.model$mod$ETA[[1]]$b)
      r.p.phenos <- NULL
      r.p.eff <- real.pred$output.model$mod$ETA[[1]]$b
    }
    if(est_h == T){
      h <- real.pred$h
    }
    if(ABC_scheme != "C"){
      rm(real.pred)
      r.p.phenos <- NULL
      r.p.eff <- NULL
    }
  }
  else{
    r.p.phenos <- NULL
  }

  # run the ABC
  ## serial
  if(par == F){

    # initialize effects storage
    out.effects <- data.table::as.data.table(matrix(NA, nrow = nrow(x), ncol = iters))
    
    for(i in 1:iters){
      cat("Iter: ", i, ".\n")
      if(is.numeric(run_number)){rn <- run_number}
      else{rn <- i}
      tout <- loop_func(x, phenotypes, out[i,"pi"], out[i,"df"], out[i,"scale"], method, ABC_scheme, t_iter = rn, r.p.phenos = r.p.phenos, r.p.eff = r.p.eff)
      out[i,"dist"] <- tout$dist
      data.table::set(out.effects, j = i,  value = tout$e)
    }
  }
  # parallel
  else{
    parms <- out[,-ncol(out)]
    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)
    
    # divide up into ncore chunks
    chunks <- split(as.data.frame(out), sample(1:par, nrow(out), replace=T))
    
    # prepare reporting function
    progress <- function(n) cat(sprintf("Chunk %d out of", n), par, "is complete.\n")
    opts <- list(progress=progress)
    
    
    output <- foreach::foreach(i = 1:par, .inorder = FALSE,
                               .options.snow = opts, 
                               .export = c("rbayesB", "get.pheno.vals", "e.dist.func", "pred", 
                                           "convert_2_to_1_column"), .packages = c("data.table", "inline"),
                               .noexport = "weighted.colSums") %dopar% {
                                 
                                 # remake the weighted.colSums function...
                                 src <- '
                                  Rcpp::NumericMatrix dataR(data);
                                  Rcpp::NumericVector weightsR(weights);
                                  int ncol = dataR.ncol();
                                  Rcpp::NumericVector sumR(ncol);
                                  for (int col = 0; col<ncol; col++){
                                  sumR[col] = Rcpp::sum(dataR( _, col)*weightsR);
                                  }
                                  return Rcpp::wrap(sumR);'
                                 
                                 weighted.colSums <- inline::cxxfunction(
                                   signature(data="numeric", weights="numeric"), src, plugin="Rcpp")
                                 
                                 
                                 out <- chunks[[i]]
                                 out.effects <- data.table::as.data.table(matrix(NA, nrow = nrow(x), ncol = nrow(out)))

                                 # run once per iter in this chunk
                                 for(i in 1:nrow(out)){
                                   if(is.numeric(run_number)){rn <- run_number}
                                   else{rn <- i}
                                   tout <- loop_func(x, phenotypes, out[i,"pi"], out[i,"df"], out[i,"scale"], method, ABC_scheme, t_iter = rn, r.p.phenos = r.p.phenos)
                                   out[i,"dist"] <- tout$dist
                                   data.table::set(out.effects, j = i,  value = tout$e)
                                 }
                                 list(dists = out, effects = out.effects)
                               }
    
    
    # release cores and clean up
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
    gc();gc()

    # bind, relies on rvest::pluck to grab only the first or only the second part
    out <- dplyr::bind_rows(rvest::pluck(output, 1))
    out.effects <- dplyr::bind_cols(rvest::pluck(output, 2))
  }

  return(list(dists = out, effects = out.effects))
}
