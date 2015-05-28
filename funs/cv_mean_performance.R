## compute crossvalidation performance of mean prediction.

rm(list = ls())

require(plyr)
require(psybayes)
require(ROCR)
source(paste0(getwd(),'/funs/helpers.R'))

load(paste0(getwd(),'/output/csf_data_subset.RData'))

mean_crossval <- data.frame()
mean_loglik <- data.frame()

cum_ps <- rep(NA, times = nrow(dat))

for (fold in 1 : length(unique(dat$cv_group))){
  train_dat <- subset(dat, cv_group != fold)
  test_dat <- subset(dat, cv_group == fold)
  
  p <- rep(mean(train_dat$correct), times = nrow(test_dat))
  
  cum_ps[dat$cv_group == fold] <- p
  
  ll <- calc_ll(y = test_dat$correct, p = p, size = 1)
  ll_norm <- ll / nrow(test_dat)
  
  pred <- prediction(predictions = p, labels = test_dat$correct)
  perf <- performance(pred, "auc")
  auc <- as.numeric(perf@y.values)
  
  this_frame <- data.frame(fold = fold, model = "mean",
                           value = auc)
  mean_crossval <- rbind(mean_crossval, this_frame)  
  
  # each log likelihood (trials):
  ll_vec <- dbinom(test_dat$correct, size = 1, prob = p, log = TRUE)
  this_frame <- data.frame(fold = fold, model = "mean",
                           lls = ll_vec)
  mean_loglik <- rbind(mean_loglik, this_frame)  
}


# overall log likelihood:
ll <- calc_ll(y = dat$correct, p = cum_ps, size = 1)
ll_norm <- ll / nrow(dat)

pred <- prediction(predictions = cum_ps, labels = dat$correct)
perf <- performance(pred, "auc")
auc <- as.numeric(perf@y.values)

perf <- performance(pred,"tpr","fpr")
plot(perf)

this_frame <- data.frame(fold = NA, model = "mean_total",
                         value = auc)
mean_crossval <- rbind(mean_crossval, this_frame)  

mean_crossval$fold <- factor(mean_crossval$fold)
mean_loglik$fold <- factor(mean_loglik$fold)

save(mean_crossval, mean_loglik, file = paste0(getwd(), "/output/mean_model_crossval.RData"))
