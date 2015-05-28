# functions for doing cross validation of GLM and transducer models. 
# Saves the matrix of p_hat for each trial.

do_crossval <- function(model, debug = FALSE){
  
  require(psybayes)
  require(rstan)
  require(plyr)
  source(paste0(getwd(),'/funs/helpers.R'))
  
  iter <- 50000
  n_saved_samples <- 1000
  
  load(paste0(getwd(),'/output/csf_data_subset.RData'))
  
  # faster debugging of fits:
  if (debug == TRUE){
    iter <- 100
    n_saved_samples <- 100
  }
  
  if (model == 1){
    model_file <- paste0(getwd(),'/funs/single_level_glm.stan')
    output_file <- paste0(getwd(),'/output/crossval_single_level_glm.RData')
  }
  
  if (model == 2){
    model_file <- paste0(getwd(),'/funs/multilevel_glm.stan')
    output_file <- paste0(getwd(),'/output/crossval_multilevel_glm.RData')
  }
  
  if (model == 3){
    model_file <- paste0(getwd(),'/funs/single_level_transducer.stan')
    output_file <- paste0(getwd(),'/output/crossval_single_level_transducer.RData')
  }
  
  if (model == 4){
    model_file <- paste0(getwd(),'/funs/single_level_transducer_B.stan')
    output_file <- paste0(getwd(),'/output/crossval_single_level_transducer_B.RData')
  }
  
  phat_matrix <- matrix(rep(NA, times = n_saved_samples*NROW(dat)), ncol = n_saved_samples)
  
  # do each fold serially (mcmc sampling will be done in parallel):
  for (fold in 1 : length(unique(dat$cv_group))){
    train_dat <- subset(dat, cv_group != fold)
    test_dat <- subset(dat, cv_group == fold)
    
    if(model == 1 | model == 2){
      # if GLM models...
      
      out_list <- make_design_matrix(train_dat, return_scale_params = TRUE)
      
      X <- out_list$X
      scale_factors <- list(center = out_list$center, scale = out_list$scale)
      
      stan_dat <- list(N = nrow(train_dat), 
                       y = train_dat$correct,
                       X = X,
                       D = ncol(X),
                       ss = as.numeric(train_dat$subject),
                       S = length(levels(train_dat$subject))) 
      
      # Load and fit stan model -------------------------
      fit <- stan_sample(file_path=model_file, stan_dat, iter = iter, n_saved_samples = n_saved_samples)
      
      # do model predictions (using training data scale factors):
      phat <- phat_samples_glm(dat = test_dat, param_list = extract(fit),
                               scale_factors = scale_factors)
    }
    
    if(model == 3 | model == 4){
      # if transducer models...
      stan_dat <- list(N = nrow(train_dat), 
                       y = train_dat$correct,
                       ped = train_dat$c_centre_target,
                       ped_plus_inc = train_dat$c_centre_target + train_dat$increment,
                       ss = as.numeric(train_dat$subject),
                       S = length(levels(train_dat$subject))) 
      
      # Load and fit stan model -------------------------
      fit <- stan_sample(file_path=model_file, stan_dat, iter = iter, n_saved_samples = n_saved_samples)
      
      # do model predictions:
      phat <- phat_samples_transducer(dat = test_dat, param_list = extract(fit))
    }
    
    # save back into phat_matrix:
    phat_matrix[which(dat$cv_group==fold),] <- phat
    save(phat_matrix, file=output_file) # save after every fold in case something dies...
  }
  
}
