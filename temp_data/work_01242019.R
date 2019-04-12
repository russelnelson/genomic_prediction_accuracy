phenos <- get.pheno.vals(x, meta$effect, h)


init <- init_pop(x = list(x = x, meta = meta, h = h, phenotypes = list(p = phenos$p)), 
                 init_gens = 3, 
                 growth.function = l_g_func, 
                 survival.function = survival.dist.func, 
                 rec.dist = rec.dist, 
                 var.theta = 0,
                 chr.length = chrl, 
                 fgen.pheno = TRUE,
                 print.all.freqs = F)


pred_vals <- pred(x = as.matrix(init$x),
                  h = .5,
                  meta = meta,
                  phenotypes = init$phenos,
                  chr.length = chrl,
                  prediction.program = "BGLR",
                  effect.sizes = meta$effect,
                  prediction.model = "BayesB",
                  make.ig = make.ig, 
                  sub.ig = 50000, 
                  maf.filt = maf, 
                  julia.path = julia.path, 
                  runID = runID, 
                  qtl_only = F,
                  pass.resid = pass.resid,
                  pass.var = pass.var, 
                  chain_length = 5000, 
                  burnin = 1000, 
                  thin = 100,
                  ntree = 1000000, 
                  null.tree = null.tree,
                  boot.ntree = boot.ntree, 
                  standardize = T,
                  save.meta = save.meta)





# pr.gs <- gs(x = pred_vals, 
#             pred.method = "model",
#             gens = max_gens,
#             growth.function = l_g_func, 
#             survival.function = survival.dist.func, 
#             selection.shift.function = sopt.func, 
#             rec.dist = rec.dist,
#             var.theta = var.theta, 
#             plot_during_progress = F, 
#             chr.length = chrl, 
#             print.all.freqs = F,
#             adjust_phenotypes = T,
#             intercept_adjust = F,
#             fgen.pheno = TRUE)
# 
# 
# 
# 
# 
# 
# 
# x2 <- convert_2_to_1_column(pred_vals$x)
# plot(x2%*%pred_vals$BGLR_mod$ETA[[1]]$b, pred_vals$a.eff$a)
# 
# 
# 
# 
# pred_vals <- pred(x = x,
#                   h = .5,
#                   meta = meta,
#                   phenotypes = NULL,
#                   chr.length = chrl,
#                   prediction.program = "BGLR",
#                   chain_length = 10000,
#                   burnin = 500,
#                   thin = 200,
#                   effect.sizes = meta$effect,
#                   prediction.model = "BayesB",
#                   make.ig = make.ig, 
#                   sub.ig = 50000, 
#                   maf.filt = maf, 
#                   julia.path = julia.path, 
#                   runID = runID, 
#                   qtl_only = F,
#                   pass.resid = pass.resid,
#                   pass.var = pass.var,
#                   ntree = ntree, 
#                   null.tree = null.tree,
#                   boot.ntree = boot.ntree, 
#                   standardize = T,
#                   save.meta = save.meta)







# testing the decline of genetic variance due to GP model over several gens

# bv <- pred.BV.from.model(pred_vals$output.model$mod, pred_vals$x, pred.method = "model", h = pred_vals$h, model.source = "BGLR", h.av = "fgen")
# 
# recom <- rand.mating(init$x, ncol(init$x)/2, meta, rec.dist, chrl, T, facet = "group")
# p2 <- get.pheno.vals(recom, meta$effect, h)
# recom_filt <- recom[pred_vals$kept.snps,]
# 
# bv2 <- pred.BV.from.model(pred_vals$output.model$mod, recom_filt, pred.method = "model", h = pred_vals$h, model.source = "BGLR", h.av = "fgen")
# 


# real data
test <- numeric(31)
test[1] <- var(get.pheno.vals(x, meta$effect, h)$a)
recom <- rand.mating(init$x, ncol(init$x)/2, meta, rec.dist, chrl, T, facet = "group")
test[2] <- var(get.pheno.vals(recom, meta$effect, h)$a)
for(i in 3:length(test)){
  print(i)
  recom <- rand.mating(recom, ncol(init$x)/2, meta, rec.dist, chrl, T, facet = "group")
  gc();gc();gc()
  p2 <- get.pheno.vals(recom, meta$effect, h)
  test[i] <- var(p2$a)
}


# bv prediction and subset
test2 <- matrix(0, ncol = 3, nrow = 31)
test2[1,1] <- var(get.pheno.vals(pred_vals$x, pred_vals$meta$effect, h)$a)
test2[1,2] <- var(pred.BV.from.model(pred_vals$output.model$mod, pred_vals$x, "model", "BGLR", pred_vals$h, h.av = "fgen")$a)
test2[,3] <- 1:31
recom2 <- rand.mating(pred_vals$x, ncol(init$x)/2, pred_vals$meta, rec.dist, chrl, T, facet = "group")
test2[2,1] <- var(get.pheno.vals(recom2, pred_vals$meta$effect, h)$a)
test2[2,2] <- var(pred.BV.from.model(pred_vals$output.model$mod, recom2, "model", "BGLR", pred_vals$h, h.av = "fgen")$a)

for(i in 3:nrow(test2)){
  print(i)
  recom2 <- rand.mating(recom2, ncol(init$x)/2, pred_vals$meta, rec.dist, chrl, T, facet = "group")
  test2[i,1] <- var(get.pheno.vals(recom2, pred_vals$meta$effect, h)$a)
  test2[i,2] <- var(pred.BV.from.model(pred_vals$output.model$mod, recom2, "model", "BGLR", pred_vals$h, h.av = "fgen")$a)
  gc();gc();gc()
}

colnames(test2) <- c("effects", "model", "gen")
test2 <- as.data.frame(test2, stringsAsFactors = F)

test2m <- reshape2::melt(test2, id.vars = "gen")
ggplot(test2m, aes(x = gen, y = value, color = variable)) + geom_point()
