## compute crossvalidation performance of a glm model fit to the log of contrast and increment.

rm(list = ls())

require(plyr)
require(psyphy)
require(psybayes)
require(ROCR)
source(paste0(getwd(),'/funs/helpers.R'))

load(paste0(getwd(),'/output/csf_data_subset.RData'))

ml_crossval <- data.frame()
ml_loglik <- data.frame()

cum_ps <- rep(NA, times = nrow(dat))

# correct for zeros by making smallest non-zero value:
dat$c_centre_target[dat$c_centre_target == 0] <- min(dat$c_centre_target[dat$c_centre_target != 0])
dat$increment[dat$increment == 0] <- min(dat$increment[dat$increment != 0])

for (fold in 1 : length(unique(dat$cv_group))){
  train_dat <- subset(dat, cv_group != fold)
  test_dat <- subset(dat, cv_group == fold)
  
  # fit a separate model to each subject:
  for (subj in 1:length(levels(dat$subject))){
    this_subj <- levels(dat$subject)[subj]
    
    this_train <- subset(dat, cv_group != fold & subject == this_subj)
    this_test <- subset(dat, cv_group == fold & subject == this_subj)
    
    # since the mafc.logit link optimiser is sensitive to starting values, initialise from normal glm parameter values:
    fit_1 <- glm(correct ~ scale(log(c_centre_target)) * scale(log(increment)), data = this_train, 
               family = binomial(), maxit = 1000)    
    
    fit <- glm(correct ~ scale(log(c_centre_target)) * scale(log(increment)), data = this_train, 
               family = binomial(mafc.logit(4)),
               start = coefficients(fit_1), maxit = 1000)
    
    p <- predict(fit, newdata = this_test, type = "response")

    cum_ps[dat$cv_group == fold & dat$subject == this_subj] <- p
  }
  
  # do calc for each fold:
  fold_p <- cum_ps[dat$cv_group == fold]
  ll <- calc_ll(y = test_dat$correct, p = fold_p, size = 1)
  ll_norm <- ll / nrow(test_dat)
  
  pred <- prediction(predictions = fold_p, labels = test_dat$correct)
  perf <- performance(pred, "auc")
  auc <- as.numeric(perf@y.values)
  
  this_frame <- data.frame(fold = fold, model = "ML",
                           value = auc)
  ml_crossval <- rbind(ml_crossval, this_frame)  
  
  # each log likelihood (trials):
  ll_vec <- dbinom(test_dat$correct, size = 1, prob = fold_p, log = TRUE)
  this_frame <- data.frame(fold = fold, model = "ML",
                           lls = ll_vec)
  ml_loglik <- rbind(ml_loglik, this_frame)  
}


# overall log likelihood:
ll <- calc_ll(y = dat$correct, p = cum_ps, size = 1)
ll_norm <- ll / nrow(dat)

pred <- prediction(predictions = cum_ps, labels = dat$correct)
perf <- performance(pred, "auc")
auc <- as.numeric(perf@y.values)

perf <- performance(pred,"tpr","fpr")
plot(perf)

this_frame <- data.frame(fold = NA, model = "ML_total",
                         value = auc)
ml_crossval <- rbind(ml_crossval, this_frame)  

ml_crossval$fold <- factor(ml_crossval$fold)
ml_loglik$fold <- factor(ml_loglik$fold)

ml_log_crossval <- ml_crossval
ml_log_loglik <- ml_loglik

# save(ml_crossval, ml_loglik, file = paste0(getwd(), "/output/ML_model_crossval.RData"))


# compare non-crossval:
fit_1 <- glm(correct ~ scale(log(c_centre_target)) * scale(log(increment)), data = subset(dat, subject == "S1"), 
             family = binomial(), maxit = 1000)    

log_fit <- glm(correct ~ scale(log(c_centre_target)) * scale(log(increment)), data = subset(dat, subject == "S1"), 
           family = binomial(mafc.logit(4)),
           start = coefficients(fit_1), maxit = 1000)

# compare AIC to non-log model:
fit_1 <- glm(correct ~ scale(c_centre_target) * scale(increment), data = subset(dat, subject == "S1"), 
             family = binomial(), maxit = 1000)    

lin_fit <- glm(correct ~ scale(c_centre_target) * scale(increment), data = subset(dat, subject == "S1"), 
           family = binomial(mafc.logit(4)),
           start = coefficients(fit_1), maxit = 1000)
