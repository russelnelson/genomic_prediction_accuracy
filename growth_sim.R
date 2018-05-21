#=======function to do growth and selection=======
gs <- function(x, effect.sizes, h, gens, growth.function, survival.function, 
               selection.shift.function, rec.dist,
               meta = NULL, plot_during_progress = FALSE, facet = "group", chr.length = 10000000){
  cat("Converting x to a data.table...\n")
  x <- data.table::as.data.table(x)
  
  
  #========faster colSums funciton"
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
  
  #==========function to figure out the total effect value of each individual=======
  #x: input gene copies. Individuals are assumed to be given their gene copies in order (ind one is columns one and two)
  #effect sizes: the phenotypic effect size of each loci, additive.
  #h: what the initial addative genetic variance is multiplied to give the enivronmental effect (heritability)
  #hist.a.var: numeric or "fgen", default "fgen". If numeric, gives the adative genetic variance of the initial pop in the first gen. If fgen, this assumes that this is the first gen.
  get.pheno.vals <- function(x, effect.sizes, h, hist.a.var = "fgen"){
    
    #function to add an enivronmental effect.
    e.dist.func <- function(A1, hist.a.var, h){
      esd <- sqrt(h*hist.a.var)
      env.vals <- rnorm(length(A1), 0, esd)
    }
    
    
    
    
    #get effect of each individual:
    a <- weighted.colSums(as.matrix(x), effect.sizes)

    a.ind <- a[seq(1, length(a), by = 2)] + a[seq(2, length(a), by = 2)] #add across both gene copies.
    
    #add environmental variance
    if(hist.a.var == "fgen"){
      pheno <- a.ind + e.dist.func(a.ind, var(a.ind), h)
      return(list(p = pheno, a = a.ind))
    }
    else{
      pheno <- a.ind + e.dist.func(a.ind, hist.a.var, h)
      return(list(pheno = pheno, a = a.ind))
    }
  }

  #===========prepare the first gen=========
  
  cat("Getting starting addative genetic and phenotypic values:")
  pheno <- get.pheno.vals(x, effect.sizes, h)
  
  a <- pheno$a #addative genetic values
  pheno <- pheno$p #phenotypic values
  
  cat("\n\n===============done===============\n\nStarting parms:\n\tmean genetic value:", mean(a), 
      "\n\tmean phenotypic value:", mean(pheno), "\n\taddative genetic variance:", var(a), "\n\th:", h)
  
  
  #make output matrix and get initial conditions
  out <- matrix(NA, nrow = gens + 1, ncol = 5)
  colnames(out) <- c("N", "mu_pheno", "mu_a", "opt", "diff")
  N <- ncol(x)/2 #initial pop size
  h.av <- var(a) #get the historic addative genetic variance.
  h.pv <- var(pheno) #historic phenotypic variance.
  opt <- mean(a) #starting optimal phenotype
  
  out[1,] <- c(N, mean(pheno), mean(a), opt, 0) #add this and the mean initial additive genetic variance
  
  if(plot_during_progress){
    pdat <- melt(out)
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
  
  
  #===========loop through each additional gen, doing selection, survival, and fisher sampling of survivors====
  
  
  cat("\nBeginning run...\n\n================================\n\n")
  for(i in 2:(gens+1)){
    #=========survival====
    #survival:
    s <- rbinom(out[(i-1),1], 1, survival.function(pheno, opt, h.pv)) #survive or not? Number of draws is the pop size in prev gen, surival probabilities are determined by the phenotypic variance and optimal phenotypein this gen.
    
    #if the population has died out, stop.
    if(sum(s) <= 1){
      out[i,1] <- 0
      out[i,2:4] <- rep(NA, 3)
      return(out[1:i,])
    }
    
    #what is the pop size after growth?
    out[i,1] <- round(growth.function(sum(s)))
    
    #make a new x with the survivors
    x <- x[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors
    
    #=========get parents and assign gcs for the next gen====
    #make a new x with individuals in next gen
    ##find parents
    mates <- matrix(sample(ncol(x)/2, out[i,1]*2, T), ncol = 2) #p1 and p2 are columns
    selfings <- which(mates[,1] == mates[,2]) #any selfing?
    while(length(selfings) > 0){ #correct selfing
      mates[selfings,] <- sample(ncol(x)/2, length(selfings)*2, T) #get new parents
      selfings <- which(mates[,1] == mates[,2]) #any selfing remaining?
    }
    
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
    r <- rmultinom(nrow(mates)*2*length(uf), size = 1, rec.dist) #for each passed chr (one for each chr for each parent), how many recombination events?
    num.rec <- which(r == 1) %% length(rec.dist) #get this in vector form. The number is the number of recombination events per passed chr.
    num.rec[num.rec == 0] <- length(rec.dist) # fix zeros, which should be the maximum number possible
    num.rec <- num.rec - 1 #since the first entry is no recombination
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
    
    rm(x, c1.mat, c2.mat)
    
    #=============adjust selection, get new phenotype scores, get ready for next gen====
    x <- x.next
    
    #get phenotypic/genetic values
    pheno <- get.pheno.vals(x, effect.sizes, h, hist.a.var = h.av)
    a <- pheno$a #addative genetic values
    pheno <- pheno$p #phenotypic values
    
    #adjust selection optima
    opt <- selection.shift.function(opt, iv = sqrt(h.av))
    
    #save
    out[i,2] <- mean(pheno)
    out[i,3] <- mean(a)
    out[i,4] <- opt
    out[i,5] <- opt - mean(a)
    cat("gen:", i-1, "\topt_s:", out[i-1,4], "\tend mean a:", out[i,3], "\tnum_surv:", sum(s), "\tpop size:", out[i,1],"\n")
    
    if(plot_during_progress){
      pdat <- melt(out)
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
    gc()
  }
  
  return(out)
}
