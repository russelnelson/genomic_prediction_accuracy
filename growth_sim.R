#======internal functions=======
#get phenotypic values given genotypes, effect sizes, and heritabilities. If hist.a.var is true, uses the amount of genomic variability this gen and h to figure out how big of an env effect to add. Otherwise uses the provided value (probably that in the first generation).
get.pheno.vals <- function(x, effect.sizes, h, hist.a.var = "fgen", standardize = FALSE){
  #function to add an enivronmental effect.
  e.dist.func <- function(A1, hist.a.var, h){
    esd <- sqrt((hist.a.var/h)-hist.a.var) # re-arrangement of var(pheno) = var(G) + var(E) and h2 = var(G)/var(pheno)
    env.vals <- rnorm(length(A1), 0, esd)
    return(env.vals)
  }
  

  
  #get effect of each individual:
  a <- weighted.colSums(as.matrix(x), effect.sizes) # faster than t(x)%*%effect.sizes!
  
  a.ind <- a[seq(1, length(a), by = 2)] + a[seq(2, length(a), by = 2)] #add across both gene copies.
  
  #standardize the genetic variance if requested.
  if(standardize){
    a.ind <- a.ind/sd(a.ind)
  }
  
  #add environmental variance
  if(hist.a.var == "fgen"){
    pheno <- a.ind + e.dist.func(a.ind, var(a.ind), h)
    return(list(pheno = pheno, a = a.ind))
  }
  else{
    pheno <- a.ind + e.dist.func(a.ind, hist.a.var, h)
    return(list(pheno = pheno, a = a.ind))
  }
}

#converts 2 column to 1 column genotypes and transposes
convert_2_to_1_column <- function(x){
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
gs <- function(x, 
               effect.sizes = NULL, 
               h,
               gens, 
               growth.function, 
               survival.function, 
               selection.shift.function, 
               rec.dist,
               meta,
               var.theta = 0,
               method = "effects",
               model = "observed",
               h_est = NULL,
               pred.mod = NULL,
               plot_during_progress = FALSE, 
               facet = "group", chr.length = 10000000,
               fgen.pheno = FALSE,
               intercept_adjust = FALSE,
               print.all.freqs = FALSE,
               adjust_phenotypes = FALSE,
               do.sexes = TRUE){
  cat("Initializing...\n")

  #=================checks========
  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }
  
  if(!is.null(effect.sizes)){
    if(nrow(x) != length(effect.sizes) | nrow(x) != nrow(meta)){
      stop("Provided x, effect sizes, and meta must all be of equal length!")
    }
  }
  else{
    if(nrow(x) != nrow(meta)){
      stop("Provided x, effect sizes, and meta must all be of equal length!")
    }
  }
  
  
  
  if(!method %in% c("model", "effects")){
    stop("Method must be provided. Options:\n\tmodel: predict phenotypes directly from the model provided.\n\teffects: predict phenotypes from estimated effect sizes.")
  }
  if(method == "model"){
    # will need this!
    e.dist.func <- function(A1, hist.a.var, h){
      esd <- sqrt((hist.a.var/h)-hist.a.var) # re-arrangement of var(pheno) = var(G) + var(E) and h2 = var(G)/var(pheno)
      env.vals <- rnorm(length(A1), 0, esd)
      return(env.vals)
    }
    if(!model %in% c("JWAS", "BGLR", "RF")){
      stop("To predict from the model, a JWAS, BGLR, or RF(ranger) model must be provided.\n")
    }
  }
  else{
    if(model == "RF"){
      stop("RF does not estimate effect sizes, so prediction must be done using the RF model.\n")
    }
  }
  
  
  
  #================
  # before doing anything else, go ahead and remove any loci from those provided with no effect! Faster this way.
  if(!is.null(effect.sizes)){
    if(any(effect.sizes == 0)){
      n.eff <- which(effect.sizes == 0)
      x <- x[-n.eff,]
      meta <- meta[-n.eff,]
      effect.sizes <- effect.sizes[-n.eff]
    }
  }
  
  #===========prepare the first gen=========
  # get starting phenotypes and addative genetic values
  if(length(fgen.pheno) != ncol(x)/2){
    cat("Generating starting phenotypic values from data.")
    pheno <- get.pheno.vals(x, effect.sizes, h)
    
    a <- pheno$a #addative genetic values
    pheno <- pheno$pheno #phenotypic values
  }
  else{
    cat("Using provided phenotypic values.")
    pheno <- fgen.pheno #provded phenotypic values.
    
    #if we are given effect sizes, grab a from those
    if(!is.null(effect.sizes)){
      a <- get.pheno.vals(x, effect.sizes, h)$a # genetic values from data
    }
    #otherwise pull predicted values if given
    #WORKING HERE
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #get the adjustment factor to use based on provided vs. predicted phenotypes.
    #to do this, calculate adjust_phenotypes phenotypic variances created by sampling enivironmental effects.
    #then, divide the observed phenotypic var by the mean of these.
    #to adjust, multiply future phenotypes by the square root of these values, then adjust the mean back to the correct mean.
    if(adjust_phenotypes != FALSE){
      ad.factor <- var(a)/h # expected phenotypic variance given the allelic variance.
      pred.hist.p.var <- ad.factor # save this for later
      ad.factor <- var(pheno)/ad.factor # how much do we need to adjust?
    }
    
    #if requested, get the amount to adjust phenotypes by in future gens.
    if(intercept_adjust){
      i.adj <- mean(pheno)
    }
  }
  
  #starting optimal phenotype, which is the starting mean phenotypic value.
  opt <- mean(pheno) #optimum phenotype
  
  cat("\n\n===============done===============\n\nStarting parms:\n\tstarting optimum phenotype:", opt, 
      "\n\tmean phenotypic value:", mean(pheno), "\n\taddative genetic variance:", var(a), "\n\tphenotypic variance:", var(pheno), "\n\th:", h, "\n")
  
  #make output matrix and get initial conditions
  out <- matrix(NA, nrow = gens + 1, ncol = 7)
  colnames(out) <- c("N", "mu_pheno", "mu_a", "opt", "diff", "var_a", "stochastic_opt")
  N <- ncol(x)/2 #initial pop size
  h.av <- var(a) #get the historic addative genetic variance.
  h.pv <- var(pheno) #historic phenotypic variance.

  out[1,] <- c(N, mean(pheno), mean(a), opt, 0, h.av, opt) #add this and the mean initial additive genetic variance
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
  
  #===========loop through each additional gen, doing selection, survival, and fisher sampling of survivors====
  
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
      return(out)
    }
    
    #what is the pop size after growth?
    out[i,1] <- round(growth.function(sum(s)))
    
    #make a new x with the survivors
    x <- x[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors
    
    # # check phenotypic variance...
    # temp <- get.pheno.vals(x, effect.sizes, h, hist.a.var = h.av)
    # ptemp <- data.frame(val = c(a, temp$a), class = c(rep("T0", length(a)), rep("T1", length(temp$a))))
    # temp <- temp$pheno
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
    x <- rand.mating(x, out[i,1], meta, rec.dist, chr.length, do.sexes, facet)
    
    #get phenotypic/genetic values
    if(method == "effects"){
      pheno <- get.pheno.vals(x, effect.sizes, h, hist.a.var = h.av)
      a <- pheno$a #addative genetic values
      pheno <- pheno$pheno #phenotypic values
    }
    else{
      if(model == "RF"){
        x <- as.matrix(x)
        x.c <- convert_2_to_1_column(x)
        colnames(xc) <- pred.mod$forest$independent.variable.names
        a <- predict(pred_vals$rj, as.data.frame(x2))
        pheno <- a + e.dist.func(a, h.av, h_est)
      }
    }
    
    
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

  if(print.all.freqs){
    a.fqs <- cbind(meta, a.fqs, stringsAsFactors = F)
    out <- list(summary = out, frequencies = a.fqs)
  }
  
  return(out)
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



#=======function to take in input data, use genomic prediction to estimate effect sizes based on phenotypes==============
#arugments:
#    pass.resid: NULL or numeric >= 0. A numeric value tells the function to pass the 
#                estimated residual variance in the model on to JWAS. A numeric value of 0 passes the exact variance,
#                a numeric value other than zero will fudge the variance number by up to the proportion given (1 fudges up to 100%).
#    pass.var: NULL or numeric >= 0. Like pass.resid, but for the true genetic variance.
#    standardize: Boolean. Should the addative genetic values be centered and scaled between -1 and 1 prior to entry into JWAS? Phenotypic values still won't be centered!
pred <- function(x, effect.sizes, ind.effects, chr.length = 10000000,
                 method = "JWAS",
                 chain_length = 100000, 
                 burnin = 50000,
                 thin = 100,
                 model = NULL,
                 make.ig = FALSE, sub.ig = FALSE, maf.filt = 0.05, 
                 julia.path = "julia", runID = "r1", qtl_only = FALSE,
                 pass.resid = NULL, pass.var = NULL, 
                 ntree = 500,
                 null.tree = 100,
                 boot.ntree = 500,
                 standardize = FALSE,
                 save.meta = TRUE){
  
  #============sanity checks================================
  # check that all of the required arguments are provided for the model we are running
  if(method %in% c("JWAS", "BGLR", "PLINK", "TASSEL", "RF")){
    
    # JWAS checks
    if(method == "JWAS"){
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
      
      # JWAS model-for now, only RR-BLUP.
      if(!model %in% c("RR-BLUP")){
        stop("Invalid JWAS model.")
      }
    }
    
    # BGLR checks:
    else if(method == "BGLR"){
      # chain length and burnin
      if(!is.numeric(c(chain_length, burnin))){
        stop("Chain_length and burnin must be numeric values.")
      }
    }
    
    # random forest
    else if(method == "RF"){
      # model
      if(!(model %in% c("RF", "RJ"))){
        stop("For random forest, the model must be RF (random forest) or RJ (random jungle).")
      }
      
      # ntree
      if(!is.numeric(ntree)){
        stop("ntree must be an integer!")
      }
      if(ntree != floor(ntree)){
        stop("ntree must be an integer!")
      }
      
      # null distribution checks
      # if doing a null:
      if(!is.null(null.tree)){
        
        # if it's not a matrix or array
        if(!(is.matrix(null.tree) | is.array(null.tree))){
          if(!is.numeric(null.tree)){
            stop("null.tree must either be provided or an integer.")
          }
          if(null.tree != floor(null.tree)){
            stop("null.tree must either be provided or an integer.")
          }
          else{
            if(!is.numeric(boot.ntrees)){
              stop("boot.ntrees must be an integer.")
            }
            if(boot.ntrees != floor(boot.ntrees)){
              stop("boot.ntrees must be an integer.")
            }
          }
        }
      }
      
    }
    
  }
  else{
    stop("Invalid method provided. Options: JWAS, BGLR, PLINK, TASSEL, or RF.")
  }
  
  cat("Preparing model inputs...\n")
  
  #============subfunctions=================================
  # do filters if requested
  filter_snps <- function(x, qtl_only, sub.ig, maf.filt, effect.sizes){
    # qtl_only
    if(qtl_only){
      if(sub.ig != FALSE | maf.filt != FALSE){
        warning("No subsampling (sub.ig) or maf filtering (maf.filt) will occur when qtl_only = TRUE.\n")
      }
      s.markers <- which(effect.sizes != 0)
      x <- x[s.markers,]
      meta <- meta[s.markers,]
    }
    
    # otherwise, other filters to check
    else{
      #filter low minor allele frequencies if requested (like many prediction studies will do!)
      if(maf.filt != FALSE){
        m.keep <- matrixStats::rowSums2(x)/ncol(x) >= maf.filt
        x <- x[m.keep,]
        meta <- meta[m.keep,]
      }
      
      #subset markers
      if(sub.ig != FALSE){
        if(nrow(x) > sub.ig){
          s.markers <- sort(sample(nrow(x), sub.ig))
          x <- x[s.markers,]
          meta <- meta[s.markers,]
        }
        else{
          warning("Fewer markers than sub.ig. Running all markers.\n")
        }
      }
    }
    return(list(x = x, meta = meta))
  }
  
  # convert two column to one column genotypes and transpose
  convert_1_to_2_column <- function(x){
    ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
    ind.genos <- matrix(ind.genos, nrow = ncol(x)/2, byrow = T) # rematrix and transpose!
    return(ind.genos)
  }
  
  
  #============prepare directories and other things=========
  r.ind.effects <- ind.effects # backup the ind effects.
  
  # create the directory to store results if it doesn't already exist and move over to it.
  if(!dir.exists(runID)){
    dir.create(runID)
  }
  owd <- getwd()
  setwd(runID)
  
  #============format data for prediction/GWAS==============
  
  if(method == "JWAS"){
    # make an individual effect file.
    ind.effects <- cbind(samp = as.character(1:500), phenotypes = ind.effects$p)
    write.table(ind.effects, "ie.txt", quote = F, col.names = T, row.names = F)
    
    # make an individual genotype file if it isn't already constructed.
    if(make.ig){
      #filter
      x <- filter_snps(x, qtl_only, sub.ig, maf.filt, effect.sizes)
      x <- x$meta
      x <- x$x
      
      # convert format
      ind.genos <- convert_1_to_2_column(x)
      ind.genos <- cbind(samp = 1:nrow(ind.genos), ind.genos) # add sample info
      colnames(ind.genos) <- c("samp", paste0("m", 1:(ncol(ind.genos)-1)))
      write.table(ind.genos, "ig.txt", quote = F, col.names = T, row.names = F)
    }
  }
  
  else if(method == "BGLR"){
    # filter
    x <- filter_snps(x, qtl_only, sub.ig, maf.filt, effect.sizes)
    meta <- x$meta
    x <- x$x
    
    # convert
    ind.genos <- convert_1_to_2_column(x) # rows are individuals, columns are SNPs
    colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
    rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS
    
    # prepare ETA
    ETA <- list(list(X = ind.genos, model = model, saveEffects = T))
  }
  
  else if(method == "RF"){
    x <- filter_snps(x, qtl_only, sub.ig, maf.filt, effect.sizes)
    meta <- x$meta
    x <- x$x
    
    t.x <- convert_1_to_2_column(x) # add sample info
    colnames(t.x) <- paste0("m", 1:ncol(t.x))
    
    t.eff <- data.frame(phenotype = ind.effects$pheno, stringsAsFactors = F)
    t.eff <- cbind(t.eff, t.x)
  }
  
  
  #=========run genomic prediction or GWAS and return results==========
  # for JWAS:
  if(method == "JWAS"){
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
    julia.call <- paste0(julia.call, " ", rv, " ", gv, " ", model)
    system(julia.call)
    
    #=========grab output and modify it to give the estimated effect size per locus=============
    e.eff <- read.table("est_effects.txt", header = F, sep = "\t")
    h <- read.table("h.txt")
    
    #save metadata for the selected markers if requested.
    if(save.meta){
      write.table(meta, "est_meta.txt", sep = "\t", quote = F, col.names = T, row.names = F)
    }
    
    setwd(owd)
    
    return(list(x = x, e.eff = e.eff, a.eff = r.ind.effects, meta = meta, h = as.numeric(h)))
  }
  
  # for BGLR:
  else if(method == "BGLR"){
    cat("Calling BGLR.\n")
    BGLR_mod <- BGLR::BGLR(y = ind.effects$pheno, ETA = ETA, nIter = chain_length, burnIn = burnin, thin = thin)
    
    # grab h2 estimate
    B <- BGLR::readBinMat('ETA_1_b.bin')
    h2 <- rep(NA,nrow(B))
    varU <- h2
    varE <- h2
    for(i in 1:length(h2)){
      u <- ind.genos%*%B[i,]	
      varU[i] <- var(u)
      varE[i] <- var(ind.effects$pheno-u)
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
    return(list(x = x, e.eff = e.eff, a.eff = r.ind.effects, meta = meta, h = h2))
  }
  
  # for randomForest
  else if(method == "RF"){
    cat("Running randomforest (ranger rj implementaiton).\n")
    # run the randomForest/jungle
    if(ncol(t.eff) - 1 >= 10000){
      rj <- ranger::ranger(dependent.variable.name = "phenotype", data = t.eff, num.trees = ntree, importance = "impurity", verbose = T, save.memory = T)
    }
    else{
      rj <- ranger::ranger(dependent.variable.name = "phenotype", data = t.eff, num.trees = ntree, importance = "impurity", verbose = T)
    }
    
    # if using a null distribution to scale importance values, do that:
    if(!is.null(null.tree)){
      
      # make a null distribution if not provided
      # note, should only need to do this once for each input dataset (Ne, effect.size, seq res, ect), since the null dists should converge! Don't need to do one for each individual trial.
      if(is.numeric(null.tree) & length(null.tree) == 1){
        
        # get bootstrapped phenotypes.
        null.dat <- matrix(sample(t.eff$phenotype, null.tree*length(t.eff$phenotype), replace = T), 
                           nrow = nrow(t.x), ncol = null.tree) # permuted phenotypes.
        
        # initialize null distribution storage and run simulations for RF and RJ
        null.mat <- matrix(0, null.tree, nrow(x)) # output. rows are simulations, columns are SNPs
        
        cat(paste0("Generating null distribution of size ", null.tree, ".\n\trep: 1\n"))
        for(i in 1:null.tree){
          if((i %% 10) == 0){cat("\t", i, ".\n")}
          tdat <- as.data.frame(cbind(phenotype = null.dat[,i], t.x), stringsAsFactors = F)
          trf <- ranger::ranger(phenotype ~ ., data = tdat, num.trees = boot.ntrees, importance = "impurity")
          null.mat[i,] <- trf$variable.importance
        }
        
        null.size <- null.tree
      }
      else{
        null.size <- ncol(null.tree)
      }
      p1 <- 1 - (rowSums(rj$variable.importance > t(null.mat))/null.size)
      rj$bootstrap.pval <- p1
    }
    return(list(x = x, rj = rj, a.eff = r.ind.effects, meta = meta))
  }
}

