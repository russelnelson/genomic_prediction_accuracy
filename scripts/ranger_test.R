library(ranger); library(ggplot2);
dat <- read.table("../genomic_prediction_accuracy/ML/cat_stats_ABC_.999pi_.75h_5df_scale1_scale_not_fixed_sim_gen.gz")
colnames(dat)[1:4] <- c("pi", "df", "scale", "h")
dat$pi <- log10(1 - dat$pi)
dat.o <- dat
dat <- na.omit(dat.o)

test_samps <- sample(nrow(dat), nrow(dat)*.25, replace = F)

dat_test <- dat[test_samps,]
dat_train <- dat[-test_samps,]


# ranger, works quite well even with 1000 trees only
rf1 <- ranger::ranger(data = dat_train[,-c(2:4)],
                      num.trees = 1000, mtry = ncol(dat) - 4,
                      dependent.variable.name = "pi", num.threads = 4, keep.inbag = TRUE)

## prediction error estimation via the cross validation
pe <- forestError::quantForestError(forest = rf1, X.train = dat_train[,-c(1:4)],
                                    X.test = dat_test[,-c(1:4)], Y.train = dat_train$pi, n.cores = 4)
pe$estimates$real <- dat_test$pi
pe$estimates$in_error <- ifelse(pe$estimates$real <= pe$estimates$upper_0.05 & pe$estimates$real >= pe$estimates$lower_0.05, 1, 0)

### residuals across parameter space
pe$estimates$scale <- dat_test$scale
pe$estimates$pi <- dat_test$pi
pe$estimates$df <- dat_test$df
pe$estimates$err <- pe$estimates$real - pe$estimates$pred
pe$estimates$squared.err <- pe$estimates$err^2
summary(lm(squared.err ~ df + scale*real, pe$estimates)) # note that scale  doesn't influence prediction error after accounting for the relationship between scale and pi. DF does, though: higher DFs have smaller prediction error




temp <- rf1$predictions$qerror(c(.025, 1 - .025), 1)
temp2 <- rf1$predictions$qerror(seq(0 + .001, 1 - .001, by = .001), 1)
ed <- data.frame(q = as.numeric(temp2), val = seq(0 + .001, 1 - .001, by = .001))
ed$q <- ed$q + rf1$predictions$estimates$pred[1]
ed$val[ed$val > .5] <- 1 - ed$val[ed$val > .5]
plot(rf1$predictions$estimates$pred[1] - ed$q, ed$val)

rf1$predictions$qerror(c(.05, .95), 1)


### diagnostic plots
ggplot(pe$estimates, aes(x = pred, y = real)) + geom_point(aes(color = scale)) +
  theme_bw() + scale_color_viridis_c() + geom_ribbon(aes(ymin = lower_0.05, ymax = upper_0.05, x = pred), fill = "grey", alpha = 0.5)
ggplot(pe$estimates, aes(x = real, y = in_error)) + geom_smooth(method = "glm", method.args = list(family = "binomial"))
ggplot(pe$estimates, aes(x = df, y = err)) + geom_smooth() # low dfs are correlated with over-estimations in pi (less effect sites)


saveRDS(list(forest = rf1, train = dat_train, test = dat_test, predictions = pe), "../genomic_prediction_accuracy/ML/pi_random_forest_test.RDS")



# df rf test -- terrible
rf1 <- ranger::ranger(data = dat_train[,-4],
                      num.trees = 1000, mtry = ncol(dat) - 2,
                      dependent.variable.name = "df", num.threads = 4, keep.inbag = TRUE)
pe <- forestError::quantForestError(forest = rf1, X.train = dat_train[,-c(2,4)],
                                    X.test = dat_test[,-c(2,4)], Y.train = dat_train$df, n.cores = 4)

pe$estimates$real <- dat_test$df
pe$estimates$pi <- dat_test$pi
ggplot(pe$estimates, aes(x = pred, y = real)) + geom_point(aes(color = pi)) +
  theme_bw() + scale_color_viridis_c() + geom_ribbon(aes(ymin = lower_0.05, ymax = upper_0.05, x = pred), fill = "grey", alpha = 0.5)
