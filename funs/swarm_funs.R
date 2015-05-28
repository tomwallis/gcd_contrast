# function to fit particle swarm optimiser.

swarm_fit <- function(data_set, model, parallel = TRUE, debug = FALSE, fit_training_set = TRUE){
  # load appropriate data set:
  if (data_set == 'real'){
    load(paste0(getwd(),'/output/csf_data_reduced.RData'))
    output_file <- paste0(getwd(),'/output/swarm_out_',data_set,'_model_',model,'.RData')
    
  } else {
    load(paste0(getwd(),'/output/simulated_data_',data_set,'.RData'))
    # change dat$correct to sim_correct (don't have to dig through mle_funs and put conditional flags):
    dat$correct <- dat$sim_correct
    
    output_file <- paste0(getwd(),'/output/swarm_out_sim_',data_set,'_model_',model,'.RData')
  }
  
  library(hydroPSO)
  source(paste0(getwd(),'/funs/mle_funs.R'))
  
  if (debug == TRUE){
    parallel <- 'none'
    n_cores <- 1
    maxit <- 50
    subjs_end <- 1    
  } else {
    n_cores <- 5
    maxit <- 5000
    parallel <- 'parallel'
    subjs_end <- length(levels(dat$subject))
  }
  
  if (model == 1){
    param_ranges <- list(rmax = c(0.001, 30),
                         p = c(0.2, 6),
                         q = c(0, 2),
                         z = c(0.0001, 1))
    
    # set up lower and upper bounds for params:
    lower <- c(param_ranges$rmax[1],
               param_ranges$p[1], 
               param_ranges$q[1],
               rep(param_ranges$z[1], times=6))
    upper <- c(param_ranges$rmax[2],
               param_ranges$p[2], 
               param_ranges$q[2],
               rep(param_ranges$z[2], times=6))    
  }
  
  if (model == 2 | model == 3){
    param_ranges <- list(rmax = c(0.001, 30),
                         p = c(0.2, 6),
                         q = c(0, 2),
                         z = c(0.0001, 1),
                         w = c(0, 1))
      
    # set up lower and upper bounds for params:
    lower <- c(param_ranges$rmax[1],
               param_ranges$p[1], 
               param_ranges$q[1],
               rep(param_ranges$z[1], times=6),
               rep(param_ranges$w[1], times=6))
    upper <- c(param_ranges$rmax[2],
               param_ranges$p[2], 
               param_ranges$q[2],
               rep(param_ranges$z[2], times=6),
               rep(param_ranges$w[2], times=6))    
  }
  
  
  
  # Parse data into test / training set if specified ------------------
  if (fit_training_set == TRUE){
    # only fit to training set:
    dat <- subset(dat, training_set == 'Training')
  } else {
    
  }
  
  
  # fit each subject separately -----------------------
  
  fit_frame <- data.frame()
  
  for (subj in 1 : subjs_end){
    print(this_subject <- levels(dat$subject)[subj])
    sub_dat <- subset(dat, subject == this_subject)
    # remove extra levels of subject factor:
    sub_dat$subject <- factor(sub_dat$subject)
    
    # fit this subject's parameters:
    out <- hydroPSO(fn = log_lik, dat = sub_dat, model = model,
                    lower = lower, upper = upper,
                    control = list(maxit = maxit, MinMax = 'max', reltol = 1e-9, 
                                   write2disk=FALSE, normalise = TRUE,
                                   parallel = parallel, par.nnodes = n_cores))    

    # process fitted params to list:
    fitted_pars <- mat_to_param(out$par, model = model)
    
    this_frame <- expand.grid(subject = this_subject,
                              convergence = out$convergence,
                              message = out$message,
                              value = out$value,
                              rmax = fitted_pars$rmax,
                              p = fitted_pars$p,
                              q = fitted_pars$q,
                              sf = sort(unique(dat$sf)))    
    
    this_frame$sf_factor <- levels(dat$sf_factor)
    this_frame$z <- as.numeric(fitted_pars$z)
    
    if (model != 1){
      this_frame$w <- as.numeric(fitted_pars$w)  
    }
    
    fit_frame <- rbind(fit_frame, this_frame)
  }
  
  fit_frame$model <- model
  
  save(fit_frame, file = output_file)
}
