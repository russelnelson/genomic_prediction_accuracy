library(igraph)

var_names <- c("1,000", "10,000", "Ne", 
               "genome", 
               "1/1,000", "1/10,000", "major effect", "effect size distribution",
               ".5", "1", "heritability",
               "effect sizes",
               "environmental effects",
               "addative genetic effects",
               "phenotypes",
               "all SNPs", "10,000 SNPs", "100,000 SNPs", "sequencing resolution",
               "GWAS", "GP",
               "effect size estimates",
               "demographic simulations on real effect sizes",
               "demographic simulations on estimated effect sizes")

interactions <- list(length = length(var_names))

Ne1000 <- interactions[[1]] <- c("Ne")
Ne10000 <- interactions[[2]] <- c("Ne")
Ne <- interactions[[3]] <- c("genome")
genome <- interactions[[4]] <- c("sequencing resolution", "demographic simulations on real effect sizes", "addative genetic effects")
ed1.1000 <- interactions[[5]] <- c("effect size distribution")
ed1.10000 <- interactions[[6]] <- c("effect size distribution")
edme <- interactions[[7]] <- c("effect size distribution")
effect.size.dist <- interactions[[8]] <- c("effect sizes")
h.5 <- interactions[[9]] <- c("heritability")
h1 <- interactions[[10]] <- c("heritability")
heritability <- interactions[[11]] <- c("environmental effects")
effect_sizes <- interactions[[12]] <- c("addative genetic effects", "demographic simulations on real effect sizes")
E <- interactions[[13]] <- c("phenotypes")
G <- interactions[[14]] <- c("phenotypes")
P <- interactions[[15]] <- c("GWAS", "GP", "demographic simulations on real effect sizes", "demographic simulations on estimated effect sizes")
all_snps <- interactions[[16]] <- c("sequencing resolution")
snps10000 <- interactions[[17]] <- c("sequencing resolution")
snps100000 <- interactions[[18]] <- c("sequencing resolution")
sequencing_res <- interactions[[19]] <- c("GWAS", "GP", "demographic simulations on estimated effect sizes")
gwas <- interactions[[20]] <- c("effect size estimates")
GP <- interactions[[21]] <- c("effect size estimates")
ef_s_est <- interactions[[22]] <- c("demographic simulations on estimated effect sizes")
ds.real <- interactions[[23]] <- c("")
ds.sim <- interactions[[24]] <- c("")

names(interactions) <- var_names

M <- matrix(0, nrow = length(var_names), ncol = length(var_names), byrow = TRUE)
colnames(M) <- var_names
rownames(M) <- var_names

for(i in 1:length(interactions)){
  M[i,colnames(M) %in% interactions[[i]]] <- 1
}

graph <- graph_from_adjacency_matrix(M)

tkplot(graph)


