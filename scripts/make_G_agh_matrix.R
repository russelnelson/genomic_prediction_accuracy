args <- commandArgs(TRUE)
dat <- as.character(args[1])
out <- as.character(args[2])

ind.genos <- data.table::fread(dat)
ind.genos <- as.matrix(ind.genos)
ind.genos <- ind.genos[,-1]

colnames(ind.genos) <- paste0("m", 1:ncol(ind.genos)) # marker names
rownames(ind.genos) <- paste0("s", 1:nrow(ind.genos)) # ind IDS
mig <- min(ind.genos)
G <- AGHmatrix::Gmatrix(ind.genos, missingValue = ifelse(mig == 0, NA, mig), method = "Yang", maf = .05)
colnames(G) <- rownames(ind.genos)
rownames(G) <- rownames(ind.genos)

write.table(G, out, sep = "\t", col.names = F, row.names = F, quote = F)
