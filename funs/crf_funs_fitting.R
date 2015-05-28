## functions to fit and predict contrast response functions for mcmc estimates derived from Stan.
# TSAW, contact tsawallis@gmail.com.

#-------------------------
# fitting function will fit the CRF to the data in data_frame
fit_crf <- function(x,y,model_name,iter=1000,n_chains=4,thin=5, ...){
  library(rstan)
  source(paste(getwd(),'/functions/helpers.R',sep="")) 
  # source data processing function.
  
  stanDat <- process_xy_to_stan(x,y,model_name=model_name) # assumes that x columns will be labelled.
  
  model <- find_model_function(model_name)
  
  # first fit bug test run:
  initial_fit <- stan(model_code = model, data = stanDat,iter = 10, chains = 4)
  
#   options(error = recover)
  
  # full fit:
  fit <- stan(fit = initial_fit, data = stanDat, iter = iter, chains = n_chains, verbose = F, thin=thin)
  print(fit)
  return(fit)
}

#-------------------------
# fitting function will fit the CRF to the data in data_frame
fit_crf_parallel <- function(x,y,model_name,iter=1000,n_chains=4,thin=5, cores = 4, ...){
  library(rstan)
  library(plyr)
  # set up for parallel sampling:
  library(doMC)
  options(cores=cores)
  registerDoMC()
  source(paste(getwd(),'/functions/helpers.R',sep="")) 
  # source data processing function.
  
  stanDat <- process_xy_to_stan(x,y,model_name = model_name) # assumes that x columns will be labelled.
  
  model <- find_model_function(model_name)
  
  # first fit bug test run:
  initial_fit <- stan(model_code = model, data = stanDat,iter = 10, chains = 4)
  
#   options(error = recover)
  
  #prep for parallel case
  fit_parallel_stan <- function(){
    fit <- stan(fit = initial_fit, data = stanDat, iter = iter, chains = 1, verbose = F,thin=thin)
    return(fit)
  }
  
  # full fit:
  parallel_fit <- foreach(i = 1:n_chains) %dopar% fit_parallel_stan()
  
  # squish fit objects together:
  fit <- sflist2stanfit(parallel_fit)
  
  print(fit)
  return(fit)
}

#-------------------------
# prediction function will predict new y values for a given stan fit and test_x.

# since models can contain different parameters and structures, this will contain some nasty ifelse stuff.

predict_crf <- function(stan_fit,test_x,model_name,...){
  # the test_x input must have labelled columns that the model expects (e.g., subject).
  
  # make sure helper functions are loaded:
  source(paste(getwd(),'/functions/helpers.R',sep="")) 
#   options(error = recover)
  
  # get appropriate prediction function:
  crf_sample <- sample_function(model_name)
  
  # extract mcmc samples from stan_fit object (don't do this in extract_params function to save time):
  params <- extract(stan_fit)
  n_samples <- length(params$lp__)
  
  pred_fun <- function(i){
    this_subj <- test_x$subject[i]
    subj <- as.numeric(this_subj)

    this_ped <- test_x$band_patch_Target[i]
    this_ped_plus_inc <- test_x$increment[i] + test_x$band_patch_Target[i]
    
    if(model_name=='csf_model' | model_name=='csf_model_2'){
      this_sf <- test_x$sf[i]
      param_list <- extract_params(params,model_name,subj_number=subj,sf_value=this_sf)
    } else {
      this_sf <- test_x$sf_factor[i]
      param_list <- extract_params(params,model_name,subj_number=subj,sf_number=as.numeric(this_sf))
    }
    
    # generate pedestal response for every mcmc sample:
    r_ped <- sapply(1:n_samples,crf_sample,c=this_ped,param_list)
    r_ped_plus_inc <- sapply(1:n_samples,crf_sample,c=this_ped_plus_inc,param_list)
    
    delta_r <- r_ped_plus_inc - r_ped

    # andrew's approximation for 4AFC:
    pc <- pnorm(delta_r / sqrt(2))^log2(4)
    return(y <- mean(pc))
  }
  
  return(predictions <- sapply(1:NROW(test_x),pred_fun))
}

