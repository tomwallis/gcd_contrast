## compute crossvalidation performance of glm model (max likelihood).

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

for (fold in 1 : length(unique(dat$cv_group))){
  train_dat <- subset(dat, cv_group != fold)
  test_dat <- subset(dat, cv_group == fold)
  
  # fit a separate model to each subject:
  for (subj in 1:length(levels(dat$subject))){
    this_subj <- levels(dat$subject)[subj]
    
    this_train <- subset(dat, cv_group != fold & subject == this_subj)
    this_test <- subset(dat, cv_group == fold & subject == this_subj)
    
    out_list <- make_design_matrix(this_train, return_scale_params = TRUE)
    
    X <- out_list$X
    scale_factors <- list(center = out_list$center, scale = out_list$scale)

    # since the mafc.logit link optimiser is sensitive to starting values, initialise from normal glm parameter values:
    fit_1 <- glm.fit(x = X, y = this_train$correct, 
               family = binomial())    
    
    fit <- glm.fit(x = X, y = this_train$correct,
               family = binomial(mafc.logit(4)),
               start = coefficients(fit_1))
    
    test_fit <- make_design_matrix(this_test, scale_factors = scale_factors)

    # do prediction manually to take scale factor into account:
    lin_pred <- test_fit %*% coef(fit)
    link_inv <- mafc.logit(4)$linkinv
    p <- link_inv(lin_pred)

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

save(ml_crossval, ml_loglik, file = paste0(getwd(), "/output/ML_model_crossval.RData"))
