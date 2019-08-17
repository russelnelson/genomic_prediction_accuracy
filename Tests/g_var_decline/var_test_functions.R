Odegard_pca_BV <- function(X, b, p.max, hist.V = FALSE, std = T){
  rs <- rowMeans(X)
  fixed <- which(rs == 0 | rs == 1)
  if(length(fixed) > 0){
    X <- X[-fixed,]
    if(!is.null(b)){
      b <- b[-fixed]
    }
    if(is.matrix(hist.V) == TRUE){
      hist.V <- hist.V[-fixed,]
    }
  }
  
  X <- convert_2_to_1_column(X)
  
  if(is.matrix(hist.V) == FALSE){
    e <- svd(X)
    V <- e$v
    percents <- cumsum(e$d)/sum(e$d)
    p <- max(which(percents <= p.max))
    V <- V[,1:p]
  }
  else{
    V <- hist.V
  }
  bT <- X%*%V
  if(is.null(b)){
    return(list(C = bT, g = NULL))
  }
  
  if(std){
    std.fac <- rowSums(abs(V))
    std.V <- V/std.fac
    shat <- t(std.V)%*%matrix(b)
  }
  else{
    shat <- t(V)%*%matrix(b)
  }
  ghat <- bT%*%shat
  # if(any(is.nan(ghat))){browser()}
  return(list(C = bT, g = ghat))
}


# makes a concatenated C matrix after Odegard eta al. 2018 using window sliding across a specific chromosome. Needs to be looped for multiple chromosomes.
make_window_cat_C <- function(X, X.meta, ws, p.max){
  tps <- X.meta$position #get the site positions, sort
  lsp <- max(tps) #get the position of the last site to use as endpoint
  ws <- 1000*ws
  c <- ws/2
  cat.C <- matrix(NA, ncol(X)/2, ncol(X/2)*ceiling(lsp/ws))
  i <- 1
  col.tracker <- 1
  while (c <= lsp){
    start <- c - ws/2 #change window start
    end <- c + ws/2 #change window end
    
    #take all the snps in the window, calculate T's theta, W's theta, and T's D
    wsnps <- which(tps <= end & tps > start) #get only the sites in the window
    if(length(wsnps) == 0){ #if no snps in window
      c <- c + ws #step along
      next #go to the next window
    }
    
    # subset the input data
    t.X <- X[wsnps,]
    out <- Odegard_pca_BV(t.X, NULL, p.max, hist.V = F, std = F)$C
    cat.C[,col.tracker:(col.tracker + ncol(out) - 1)] <- out 
    col.tracker <- col.tracker + ncol(out)
    
    c <- c + ws
    i <- i + 1
  }
  
  # remove excess columns!
  if(col.tracker < ncol(cat.C)){
    cat.C <- cat.C[,-c(col.tracker:ncol(cat.C))]
  }
  
  return(cat.C)
}

# makes a concatenated C matrix after Odegard eta al. 2018 using window sliding across a specific chromosome. Needs to be looped for multiple chromosomes.
make_window_cat_ghat <- function(X, X.meta, b, ws, p.max, std = T){
  tps <- X.meta$position #get the site positions, sort
  lsp <- max(tps) #get the position of the last site to use as endpoint
  ws <- 1000*ws
  c <- ws/2
  cat.ghat <- matrix(NA, ncol(X)/2, ceiling(lsp/ws))
  i <- 1
  while (c <= lsp){
    start <- c - ws/2 #change window start
    end <- c + ws/2 #change window end
    
    #take all the snps in the window, calculate T's theta, W's theta, and T's D
    wsnps <- which(tps <= end & tps > start) #get only the sites in the window
    if(length(wsnps) == 0){ #if no snps in window
      c <- c + ws #step along
      next #go to the next window
    }
    
    # subset the input data
    t.X <- X[wsnps,]
    t.b <- b[wsnps]
    out <- Odegard_pca_BV(t.X, t.b, p.max, hist.V = F, std = std)$g
    cat.ghat[,i] <- out 
    c <- c + ws
    i <- i + 1
  }
  
  return(matrixStats::rowSums2(cat.ghat))
}




test_var <- function(x, e, h, trials, p.max, model, rec.dist, x.meta, chrl, ws, do.sexes = T, facet = "group", hist.V = F, std = T){

  var_test <- matrix(NA, nrow = trials, ncol = 11)
  colnames(var_test) <- c("gen", "var_BV", "var_eBV", "var_pcaBV", "var_wsBV",
                          "cor_BV_eBV", "cor_BV_pcaBV", "cor_eBV_pcaBV", "cor_BV_wsBV",
                          "cor_eBV_wsBV", "cor_pcaBV_ws_BV")
  BVs <- array(NA, dim = c(trials, ncol(x)/2, 4))
  var_test[,1] <- 1:trials
  recom <- x
  
  e.eff <- model$e.eff$V2
  X.meta <- x.meta[model$kept.snps,]
  
  for(i in 1:trials){
    print(i)
    # true bv
    BVs[i,,1] <- get.pheno.vals(recom, e, h)$a
    var_test[i,2] <- var(BVs[i,,1])
    
    # keep the correct values
    maf.recom <- recom[model$kept.snps,]
    
    # estimated bv
    BVs[i,,2] <- pred.BV.from.model(model$output.model$mod, maf.recom, "model", "BGLR", model$h, h.av = "fgen")$a
    var_test[i,3] <- var(BVs[i,,2])
    var_test[i,6] <- cor(BVs[i,,1], BVs[i,,2])
    
    # pca bv
    BVs[i,,3] <- Odegard_pca_BV(maf.recom, e.eff, p.max, hist.V, std)$g
    var_test[i, 4] <- var(BVs[i,,3])
    var_test[i, 7] <- cor(BVs[i,,1], BVs[i,,3])
    var_test[i, 8] <- cor(BVs[i,,2], BVs[i,,3])
    
    # ws bv
    BVs[i,,4] <- make_window_cat_ghat(maf.recom, X.meta, e.eff, ws, p.max, std)
    var_test[i, 5] <- var(BVs[i,,4])
    var_test[i, 9] <- cor(BVs[i,,1], BVs[i,,4])
    var_test[i, 10] <- cor(BVs[i,,2], BVs[i,,4])
    var_test[i, 11] <- cor(BVs[i,,3], BVs[i,,4])
    
    # if(is.na(var_test[i,4])){browser()}
    
    
    if(i != nrow(var_test)){
      recom <- rand.mating(recom, ncol(recom)/2, x.meta, rec.dist = rec.dist,
                           chrl, do.sexes, facet = facet)
    }
    
    gc();gc();gc();gc();gc()
  }

  var_test <- as.data.frame(var_test)
  var_test$pvar_BV <- var_test$var_BV/var_test$var_BV[1]
  var_test$pvar_eBV <- var_test$var_eBV/var_test$var_eBV[1]
  var_test$pvar_pcaBV <- var_test$var_pcaBV/var_test$var_pcaBV[1]
  var_test$pvar_wsBV <- var_test$var_wsBV/var_test$var_wsBV[1]
  return(list(var_test = var_test, BVs = BVs))
}