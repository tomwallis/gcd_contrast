# script to pool and plot crossvalidation data.

rm(list = ls())
library(ggplot2)
library(reshape2)
library(psybayes)
library(ROCR)
library(dplyr)
library(grid)
library(gridExtra)

source(paste0(getwd(),'/funs/helpers.R'))

# load data frame:
load(paste0(getwd(),'/output/csf_data_subset.RData'))

# load mean model crossval:
load(paste0(getwd(), "/output/mean_model_crossval.RData"))
load(paste0(getwd(), "/output/ML_model_crossval.RData"))

compute_auc <- function(y, p){
  pred <- prediction(predictions = p, labels = y)
  perf <- performance(pred, "auc")
  return(auc <- as.numeric(perf@y.values))
}

# calc for each model, add to cumulative df -----------------------
log_liks <- data.frame()
auc_folds <- data.frame()
auc_total <- data.frame()

for (model in 1 : 4){
  if (model == 1){
    model_name <- 'GLM'
    model_file <- paste0(getwd(),'/output/crossval_single_level_glm.RData')
  }
  
  if (model == 2){
    model_name <- 'Multilevel GLM'
    model_file <- paste0(getwd(),'/output/crossval_multilevel_glm.RData')
  }
  
  if (model == 3){
    model_name <- 'Transducer A'
    model_file <- paste0(getwd(),'/output/crossval_single_level_transducer.RData')
  }
  
  if (model == 4){
    model_name <- 'Transducer B'
    model_file <- paste0(getwd(),'/output/crossval_single_level_transducer_B.RData')
  }
  
  load(model_file)
  
  cum_p <- rep(NA, times = nrow(dat))
  
  # for each fold, and each mcmc sample, generate a LL and ROC:
  for (fold in 1 : length(unique(dat$cv_group))){
    this_phat <- phat_matrix[dat$cv_group==fold,]
    this_y <- subset(dat, cv_group == fold)$correct
    
    # mean posterior p:
    p <- rowMeans(this_phat)
    cum_p[dat$cv_group==fold] <- p
    
    # log lik for each trial:
    ll_vec <- dbinom(this_y, size = 1, prob = p, log = TRUE)
    
#     # log likes as differences to baseline:
#     baseline <- mean_loglik$lls[dat$cv_group==fold]
#     ll_vec <- ll_vec - baseline
    
    this_out <- data.frame(fold = fold, model = model_name, 
                           value = ll_vec)
    log_liks <- rbind(log_liks, this_out)
    
    # area under rocs:
    auc <- compute_auc(y = this_y, p = p)
    
    this_out <- data.frame(fold = fold, model = model_name, 
                           value = auc)
    auc_folds <- rbind(auc_folds, this_out)
  }
  
  # do prediction on whole data set:
  # area under rocs:
  auc <- compute_auc(y = dat$correct, p = cum_p)
  
  this_out <- data.frame(model = model_name, 
                         value = auc)
  auc_total <- rbind(auc_total, this_out)
}

# Process and plot fold info -----------------------------

auc_folds$variable <- "Area under ROC"
log_liks$variable <- "Log likelihood"

combined <- rbind(auc_folds, log_liks)
combined$fold <- factor(combined$fold)
combined$variable <- factor(combined$variable)

# summarise average loglik of mean model and ML for each fold:
mean_model <- mean_loglik %.%
  group_by(fold) %.%
  summarise(value = mean(lls))
mean_model$variable <- "Log likelihood"

ml_model <- ml_loglik %.%
  group_by(fold) %.%
  summarise(value = mean(lls))
ml_model$variable <- "Log likelihood"

# attach AUC data from ml model:
blah <- subset(ml_crossval, !is.na(fold), select = -model)
blah$variable <- "Area under ROC"
ml_model <- rbind(ml_model, blah)

# plots -----------------------------
fig <- ggplot(combined, aes(x = model, y = value)) + 
  facet_grid(variable ~ fold, scales = "free_y") + 
  geom_hline(data = mean_model, aes(x = NULL, yintercept=value)) + 
  geom_hline(data = ml_model, aes(x = NULL, yintercept=value), linetype = 2) + 
  geom_violin(size = 0.2) + 
  plot_hdi_pointrange() +
  xlab("") + ylab("") +
  theme_minimal(base_size=11) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
  theme(axis.text.x = element_text(angle = 90))
# fig
ggsave(file = paste0(getwd(), "/figs/crossval_folds.pdf"), width = 7, height = 6)


# Process and plot total -----------------------------

# do separate plots for AUC and LL to allow easier axis scaling.
auc_ml <- subset(ml_crossval, model == "ML_total", select = c(-fold))

fig <- ggplot(auc_total, aes(x = model, y = value)) +
  geom_point() + 
  xlab("") + ylab("Area under ROC") + 
  theme_minimal(base_size=11) +
  geom_hline(data = auc_ml, aes(x = NULL, yintercept=value), linetype = 2) + 
  scale_y_continuous(breaks = seq(0.68, 0.71, length = 4), limits = c(0.68, 0.712))
fig_1 <- fig

ll_ml <- data.frame(value = mean(ml_loglik$lls))
ll_mean_model <- data.frame(value = mean(mean_loglik$lls))

ll_models <- log_liks %.%
  group_by(model) %.%
  summarise(value = mean(value))

# shorter model names:
levels(ll_models$model) <- c("GLM", "MultiGLM", "Trans A", "Trans B")

### Convert log likelihood to bits (log base 2), express relative to mean model:
ll_mean_model <- ll_mean_model / log(2)  

ll_models$value <- (ll_models$value / log(2)) - ll_mean_model$value
ll_ml <- (ll_ml / log(2)) - ll_mean_model$value

log_liks$rel <- ((log_liks$value) / log(2)) - ll_mean_model$value

# fig <- ggplot(ll_models, aes(x = model, y = value)) +
#   geom_point() + 
#   xlab("") + ylab("Log likelihood (bits / trial)") + 
#   theme_minimal(base_size=11) +
#   geom_hline(data = ll_ml, aes(x = NULL, yintercept=value), linetype = 2) +
#   scale_y_continuous()
# # fig
# fig_2 <- fig

# bootstrapped errors:
fig <- ggplot(log_liks, aes(x = model, y = rel)) +
  stat_summary(fun.data = "mean_cl_boot", B=5000) + 
  xlab("") + ylab("Log likelihood (bits / trial)") + 
  theme_minimal(base_size=11) +
  geom_hline(data = ll_ml, aes(x = NULL, yintercept=value), linetype = 2) +
  scale_y_continuous(breaks=c(0, .05, .1))
# fig
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

pdf(file=paste0(getwd(),'/figs/crossval_summaries.pdf'),width=3,height=3.5)
grid.arrange(fig_1, fig_2, ncol = 1, widths = c(0.5, 0.5), heights = c(0.5, 0.5))
dev.off()



# sign rank test?
# wilcox.test(log_liks$value[log_liks$model == "Transducer B"], mean_loglik$lls, paired = TRUE, alternative = "greater")

