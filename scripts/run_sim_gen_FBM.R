library(GeneArchEst)
#=========generate a bunch of sims============
# ABC params
## priors
df_func <- function(x) runif(x, 1, 100)


args <- commandArgs(TRUE)
genofile <- as.character(args[1])
ABCfile <- as.character(args[2])
metafile <- as.character(args[3])
phenofile <- as.character(args[4])
windowfile <- as.character(args[5])
gmmatfile <- as.character(args[6])
grmfile <- as.character(args[7])
hmean <- as.numeric(args[7])
hsd <- as.numeric(args[8])
outname <- as.character(args[9])

iters <- 1

# read in the genotypes
x <- bigstatsr::big_attach(genofile)

# read in the phenotypes
phenos <- readRDS(phenofile)
phenos <- phenos$phenos

# read in the ABC results
res <- readRDS(ABCfile)$ABC_res

# read in the windows
pass_windows <- readRDS(windowfile)

# read in the G matrix
pass_G <- data.table::fread(grmfile, header = F)
pass_G <- as.matrix(pass_G)

# read in the subset metadata
meta <- as.data.frame(data.table::fread(metafile, header = F))
colnames(meta) <- c("chr", "position")

# clean
gc(); gc();



# run
sims <- sim_gen(x = x, meta = meta, iters = iters, center = T, scheme = "gwas", 
                parameter_distributions = list(pi = "joint", scale = "joint", d.f = df_func), 
                h_dist = function(x) rnorm(x, hmean, hsd), joint_res = res, joint_acceptance = 0.005, 
                joint_res_dist = "ks", pass_windows = pass_windows, pass_G = pass_G, GMMAT_infile = gmmatfile)


data.table::fwrite(sims$stats, outname, sep = "\t", quote = F, col.names = F, row.names = F)


