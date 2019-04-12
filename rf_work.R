require(gradientForest)

convert_1_to_2_column <- function(x){
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  ind.genos <- matrix(ind.genos, nrow = ncol(x)/2, byrow = T) # rematrix and transpose!
  return(ind.genos)
}

meta$effect <- 0
sub <- sample(length(meta$effect), 1000, replace = F)
t.meta <- meta[sub,]
tr.x <- x[sub,]
t.meta$effect <- 0
t.meta$effect[sample(100, 20)] <- 1

t.eff <- get.pheno.vals(tr.x, t.meta$effect, .5)

t.x <- convert_1_to_2_column(tr.x) # add sample info
colnames(t.x) <- paste0("m", 1:ncol(t.x))

t.eff <- data.frame(phenotype = t.eff$pheno, stringsAsFactors = F)
t.eff <- cbind(t.eff, t.x)

# construct null distribution
null.size <- 100
null.mat <- array(0, dim = c(null.size, nrow(tr.x), 2)) # output. first index is sim num, second is SNP, third is the measure
null.dat <- matrix(sample(t.eff$phenotype, length(null.dat), replace = T), 
                   nrow = nrow(t.x), ncol = null.size) # permuted phenotypes.

#note, should only need to do this once for each input dataset (Ne, effect.size, seq res, ect), since the null dists should converge! Don't need to do one for each individual trial.
cat(paste0("Generating null distribution of size ", null.size, ".\n\trep: 1\n"))
for(i in 1:null.size){
  if((i %% 10) == 0){cat("\t", i, ".\n")}
  tdat <- cbind(phenotype = null.dat[,i], t.x)
  trf <- randomForest::randomForest(phenotype ~ ., data = tdat, ntree = 100, importance = T)
  null.mat[i,,] <- trf$importance
}

rf <- randomForest::randomForest(phenotype ~., data = t.eff, ntree = 500, importance = T)

p1 <- 1 - (rowSums(rf$importance[,1] > t(null.mat[,,1]))/null.size)
p2 <- 1 - (rowSums(rf$importance[,2] > t(null.mat[,,2]))/null.size)

# Does a good job at this:
plot(t.meta$effect*rowMeans(tr.x) ~ log(p1))



x2 <- rand.mating(pred_vals$x, N.next = 20, rec.dist = rec.dist, 
                  chr.length = chrl, meta = pred_vals$meta, 
                  do.sexes = F, facet = "group")
x2 <- as.matrix(x2)
x2 <- convert_1_to_2_column(x2)

colnames(x2) <- pred_vals$rj$forest$independent.variable.names

rjp <- predict(pred_vals$rj, as.data.frame(x2))

library(ranger)






