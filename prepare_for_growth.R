#=========import and prepare data=========
x <- "theta4k_1000_10_rho40k.txt"
chrl <- 10000000

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
    tg <- matrix(tg, ncol = chrls[i], nrow = gc, byrow = T) #put into a matrix
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

x <- process_ms(x, chrl)
meta <- x$meta
x <- x$x

#========add effect sizes=========
prob.effect <- 0.01
effect.dist.func <- function(n){
  x <- rnorm(n, 0, .5)
  # (x - min(x))/(max(x) - min(x))
}

#does each have an effect?
meta$effect <- rbinom(nrow(meta), 1, prob.effect)

#what are the effect sizes?
meta$effect[which(meta$effect == 1)] <- effect.dist.func(sum(meta$effect))

#=======get the initial phenotypic variance========
x <- matrix(as.numeric(x), nrow(x))

# 'fold' the spectrum by randomly reassigning 0s and 1s (so that the allele with effect size != 0 isn't necissarily the derived allele)
# ssf <- rbinom(nrow(x), 1, .5)
# xf <- x
# xf[which(as.logical(ssf)),] <- ifelse(xf[which(as.logical(ssf)),] == 1, 0, 1)
# rm(x)

#=======growth and survival functions===========
#logarithmic growth
l_g_func <- function(x, K = 500, r = 2){
  return((K*x*exp(r))/(K + x*(exp(r*1) - 1)))
}

#change to hist pheno var
#surivival probability follows a normal distribution around the optimal phenotype.
s_norm_func <- function(x, opt_pheno, hist_var = h.av){
  x <- dnorm(x, opt_pheno, sqrt(hist_var)*2) #normal dist, sd is equal to ancestral times 2.
  ((.75-0)*(x-min(x))/(max(x) - min(x))) #scaled between 0 and .5
  #(x-min(x))/(max(x)-min(x)) #scaled.
}

#increases the selection optimum by a percentage each gen
sopt_sp_func <- function(x, iv, slide = 0.05){
  x <- x + iv*slide
}

#increases the selection optimum according to logarithmic growth
sopt_stl_func <- function(x, K = out[1,4]*1.5, r = 3){
  if(x < 0){
    x <- abs(x)
    K <- -K
    return(-1*((K*x*exp(r))/(K + x*(exp(r*1) - 1))))
  }
  else{
    return((K*x*exp(r))/(K + x*(exp(r*1) - 1)))
  }
}


rec.dist <- 1/2^(0:floor(log(1000,2)))

# xf <- data.table::as.data.table(xf)
x <- data.table::as.data.table(x)

out <- gs(x, meta$effect, 1, 30, l_g_func, s_norm_func, sopt_sp_func, rec.dist, meta = meta, plot_during_progress = F)



N <- numeric(20)
N[1] <- -53.1794
for(i in 1:(length(N)-1)){
  N[i+1]<- sopt_stl_func(N[i], -53.1794*1.5, 1)
}
plot(N)
