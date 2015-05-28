# do crossvalidation for full glm.
rm(list = ls())
library(psybayes)
library(psyphy)
library(ROCR)
library(reshape2)
library(grid)
library(gridExtra)
library(rstan)

# load data set for this section:
load(paste0(getwd(), "/output/csf_data_subset_III.RData"))
source(paste0(getwd(), "/funs/helpers.R"))

# load stan fit:
load(paste0(getwd(),'/output/full_glm.RData'))
params <- extract(fit)

train_dat <- subset(dat, cv_group == "train")
test_dat <- subset(dat, cv_group == "test")

#-------------------------------------------------
# Mean model (baseline) test performance --------------------
#-------------------------------------------------

p <- rep(mean(train_dat$correct), times = nrow(test_dat))

ll <- calc_ll(y = test_dat$correct, p = p, size = 1)
ll_norm <- ll / nrow(test_dat)

# log lik for each trial:
ll_vec_mean <- dbinom(test_dat$correct, size = 1, prob = p, log = TRUE)

pred <- prediction(predictions = p, labels = test_dat$correct)
perf <- performance(pred, "auc")
auc <- as.numeric(perf@y.values)

mean_performance = data.frame(model='mean', auc=auc, ll=ll_norm)


#-------------------------------------------------
# Stan Model test performance --------------------
#-------------------------------------------------

# get design matrix and store training data normalisation:
out_list <- make_design_matrix_full(train_dat, return_scale_params = TRUE)
X <- out_list$X
scale_factors <- list(center = out_list$center, scale = out_list$scale)  

# generate predictions for each test trial from stan fits:
cum_ps <- trial_p_glm(dat = test_dat, param_list = params, full_model = TRUE,
                      scale_factors = scale_factors)

# log likelihood of test set:
ll <- calc_ll(y = test_dat$correct, p = cum_ps, size = 1)
ll_norm <- ll / nrow(test_dat)

# log lik for each trial:
ll_vec_multiglm <- dbinom(test_dat$correct, size = 1, prob = cum_ps, log = TRUE)

pred <- prediction(predictions = cum_ps, labels = test_dat$correct)
perf <- performance(pred, "auc")
auc <- as.numeric(perf@y.values)

perf <- performance(pred,"tpr","fpr")
plot(perf)

test_performance <- data.frame(model = "Multilevel GLM",
                         auc = auc,
                         ll = ll_norm)


#-------------------------------------------------
# ML test performance --------------------
#-------------------------------------------------

cum_ps <- rep(NA, times = nrow(test_dat))

# get design matrix and store training data normalisation:
out_list <- make_design_matrix_full(train_dat, return_scale_params = TRUE)
X <- out_list$X
scale_factors <- list(center = out_list$center, scale = out_list$scale)  

# test data design matrix:
X_test <- make_design_matrix_full(test_dat, scale_factors = scale_factors)

# annoying. have to write predict function because it returns a list?
my_predict <- function(fit, X){
  eta <- X %*% fit$coefficients
  # while I can't pass .m here, seems to have a lower asymptote of 0.25 nevertheless.
  response <- fit$family$linkinv(eta)
}

# fit a separate model to each subject:
for (subj in 1:length(levels(dat$subject))){
  this_subj <- levels(dat$subject)[subj]
  
  train_x <- X[train_dat$subject == this_subj,]
  train_y <- train_dat$correct[train_dat$subject == this_subj]
  
  test_x <- X_test[test_dat$subject == this_subj,]
  test_y <- test_dat$correct[test_dat$subject == this_subj]
  
  # since the mafc.logit link optimiser is sensitive to starting values, initialise from normal glm parameter values:
  f1 <- glm.fit(x = train_x, y = train_y, family = binomial())    
  
  f2 <- glm.fit(x = train_x, y = train_y,
                family = binomial(mafc.logit(4)), 
                start = coef(f1), control = list(maxit = 1000))
  
  p <- my_predict(f2, X = test_x)
  
  cum_ps[test_dat$subject == this_subj] <- p
}

# log likelihood of test set:
ll <- calc_ll(y = test_dat$correct, p = cum_ps, size = 1)
ll_norm <- ll / nrow(test_dat)

# log lik for each trial:
ll_vec_ml <- dbinom(test_dat$correct, size = 1, prob = cum_ps, log = TRUE)

pred <- prediction(predictions = cum_ps, labels = test_dat$correct)
perf <- performance(pred, "auc")
auc <- as.numeric(perf@y.values)

perf <- performance(pred,"tpr","fpr")
plot(perf)

this_frame <- data.frame(model = "Single ML",
                                  auc = auc,
                                  ll = ll_norm)

test_performance <- rbind(test_performance, this_frame)

print(test_performance)


### Convert log likelihood to bits (log base 2), express relative to mean model:
mean_performance$ll <- mean_performance$ll / log(2)  
test_performance$ll <- (test_performance$ll / log(2)) - mean_performance$ll

# do for trial-level data:
ll_vec_multiglm <- (ll_vec_multiglm / log(2)) - (ll_vec_mean / log(2))
ll_vec_ml <- (ll_vec_ml / log(2)) - (ll_vec_mean / log(2))

# ll_vec_multiglm <- (ll_vec_multiglm / log(2)) - mean_performance$ll
# ll_vec_ml <- (ll_vec_ml / log(2)) - mean_performance$ll

trial_lls <- expand.grid(model='Multilevel GLM', ll=ll_vec_multiglm)
this_lls <- expand.grid(model='Single ML', ll=ll_vec_ml)
trial_lls <- rbind(trial_lls, this_lls)


# plot -------------------

# test_performance <- melt(test_performance)
# levels(test_performance$variable) <- c("Area under ROC", "Log likelihood")

fig <- ggplot(test_performance, aes(x = model, y = auc)) + 
  geom_point() + 
  xlab("") + ylab("Area under ROC") + 
  theme_minimal(base_size=11) +
  scale_y_continuous(breaks = seq(0.7, 0.71, length = 3), limits = c(0.698, 0.712))
fig_1 <- fig

# fig <- ggplot(test_performance, aes(x = model, y = ll)) + 
#   geom_point() + 
#   xlab("") + ylab("Log likelihood (bits / trial)") + 
#   scale_y_continuous(breaks=seq(0, 0.1, length=3), limits = c(0, 0.1))
#   theme_minimal(base_size=11)
# fig_2 <- fig

fig <- ggplot(trial_lls, aes(x = model, y = ll)) + 
  stat_summary(fun.data='mean_cl_boot', B=5000) +
  xlab("") + ylab("Log likelihood (bits / trial)") + 
  scale_y_continuous(breaks=seq(0, 0.1, length=3)) + 
  coord_cartesian(ylim = c(0, .15)) + 
  theme_minimal(base_size=11)
fig_2 <- fig

base_size <- 8

fig_1 <- fig_1 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  ggtitle("A")

fig_2 <- fig_2 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) +
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  ggtitle("B")

pdf(file=paste0(getwd(),'/figs/crossval_full_model.pdf'),width=3,height=3.5)
grid.arrange(fig_1, fig_2, ncol = 1, widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()

save(test_performance, file = paste0(getwd(), "/output/full_glm_crossval.RData"))
