## functions for fitting CRF models via maximum likelihood to each observer's data.

mle_fit <- function(data_set, model, parallel = TRUE, debug = FALSE, fit_training_set = TRUE){
  # load appropriate data set:
  if (data_set == 'real'){
    load(paste0(getwd(),'/output/csf_data_reduced.RData'))
    output_file <- paste0(getwd(),'/output/mle_out_',data_set,'_model_',model,'.RData')
    
  } else {
    load(paste0(getwd(),'/output/simulated_data_',data_set,'.RData'))
    # change dat$correct to sim_correct (don't have to dig through mle_funs and put conditional flags):
    dat$correct <- dat$sim_correct
    
    output_file <- paste0(getwd(),'/output/mle_out_sim_',data_set,'_model_',model,'.RData')
  }
  
  n_inits <- 10
  n_cores <- 10
  maxit <- 2000
  method <- 'L-BFGS-B'
  
  if (debug == TRUE){
    trace <- 5
    parallel <- FALSE
    n_inits <- 2
    n_cores <- 1
    maxit <- 50
    method <- 'L-BFGS-B'
  } else {
    trace <- 0
  }
  
  param_ranges <- list(rmax = c(0.001, 30),
                       p = c(0.2, 6),
                       q = c(0, 2),
                       z = c(0.0001, 1),
                       w = c(0, 1))
  
  # Parse data into test / training set if specified ------------------
  if (fit_training_set == TRUE){
    # only fit to training set:
    dat <- subset(dat, training_set == 'Training')
  } else {
    
  }
  
  
  # fit each subject separately, seed inits in parallel ---------------
  
  fit_frame <- data.frame()
  
  for (subj in 1 : length(levels(dat$subject))){
    print(this_subject <- levels(dat$subject)[subj])
    sub_dat <- subset(dat, subject == this_subject)
    # remove extra levels of subject factor:
    sub_dat$subject <- factor(sub_dat$subject)
    
    # fit this subject's parameters:
    out <- optim_wrapper(dat = sub_dat, n_inits = n_inits, model = model,
                         param_ranges = param_ranges, n_cores = n_cores,
                         method = method, maxit = maxit, parallel = parallel,
                         trace = trace)

    # process fitted params to list:
    fitted_pars <- mat_to_param(out$par, model = model)
    
    this_frame <- expand.grid(subject = this_subject,
                              convergence = out$convergence,
                              value = out$value,
                              val_sd = out$val_sd,
                              n_NA = out$n_NA,
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
  
  fit_frame$method <- method
  fit_frame$maxit <- maxit
  fit_frame$model <- model
  
  save(fit_frame, file = output_file)
}



# for a given set of parameters and data, calculate likelihood of the data given the params:
log_lik <- function(par, dat, model, ...){
  require(psybayes)
  source(paste0(getwd(),'/funs/helpers.R'))
  
  # arrange data variables
  param_list <- mat_to_param(par, model)
  norm_pool <- subset(dat, select = c_surround_band_0_target : c_surround_band_5_target)
  
  # calculate p for each row:
  p_hat <- trial_p(dat = dat, param_list = param_list,
                   norm_pool = norm_pool, model = model)

  # check that no p values are 1 (will produce infs if response is 0)
  if (any(p_hat == 1)){
    warning('p_hats of 1 produced')
    p_hat[p_hat == 1] <- 0.9999
  }
  
  # calc log likelihood:
  val <- calc_ll(dat$correct, p_hat, size = 1)
  
  if (is.infinite(val)) {
    warning('likelihood calculation generated infinite values -- numerical under/overflow')
    val <- -sqrt(.Machine$double.xmax)
  }
  
  if (is.na(val)){
    warning('likelihood calculation generated NAs -- check')
    val <- -sqrt(.Machine$double.xmax)
  }
  
  # compute normalised likelihood (divide by number of trials):
  n_trials <- nrow(dat)
  val <- val / n_trials
  
  return(val) 
}


param_to_mat <- function(param_list, optim_inits = FALSE){
  
  mat <- NULL
  col_names <- NULL
  
  if (length(param_list$p) == 1){
    for(i in 1 : length(names(param_list))){
      if (names(param_list)[i] == 'w' | names(param_list)[i] == 'z'){
        j <- 1:6
        this_names <- paste0(names(param_list)[i], j)
        col_names <- c(col_names, this_names)
        mat <- cbind(mat, param_list[[names(param_list)[i]]])
      } else {
        col_names <- c(col_names, names(param_list)[i])
        mat <- cbind(mat, param_list[[names(param_list)[i]]])
      }
    }
    colnames(mat) <- col_names
    mat <- as.data.frame(mat)
    return(mat)    
  } else {
    for(i in 1 : length(names(param_list))){
      if (names(param_list)[i] == 'w' | names(param_list)[i] == 'z'){
        j <- 1:6
        this_names <- paste0(names(param_list)[i], j)
        col_names <- c(col_names, this_names)
        if (optim_inits == TRUE) {
          # calling from optim, and rows are initial vals.
          mat <- cbind(mat, param_list[[names(param_list)[i]]])  
        } else {
          mat <- cbind(mat, t(param_list[[names(param_list)[i]]]))  
        }
        
      } else {
        col_names <- c(col_names, names(param_list)[i])
        mat <- cbind(mat, param_list[[names(param_list)[i]]])
      }
    }
    colnames(mat) <- col_names
    mat <- as.data.frame(mat)
    return(mat)    
  }
}

mat_to_param <- function(mat, model){
  if (model != 1 & model != 2 & model != 3) stop("Don't know which model to fit!")
  # check that mat is a data frame:
  if (is.data.frame(mat) == FALSE){
    mat <- as.data.frame(t(mat)) # calling from optim, returns named character vector of wrong orientation.
  }

  if (model == 1) {
    
    if (names(mat)[1]=="V1" | names(mat)[1]=="Param1"){
      # calling from swarm optimiser.
      n_sfs <- 6
      names(mat) <- c('rmax','p','q', paste0('z',1:n_sfs))
    }
    
    param_list <- list(rmax = mat$rmax,
                       p = mat$p,
                       q = mat$q,
                       z = subset(mat, select = z1 : z6))
  }
  
  if (model == 2 | model == 3){
    if (names(mat)[1]=="V1" | names(mat)[1]=="Param1"){
      # calling from swarm optimiser.
      n_sfs <- 6
      names(mat) <- c('rmax','p','q', 
                      paste0('z',1:n_sfs),
                      paste0('w',1:n_sfs))
    }
    
    param_list <- list(rmax = mat$rmax,
                       p = mat$p,
                       q = mat$q,
                       z = subset(mat, select = z1 : z6),
                       w = subset(mat, select = w1 : w6))
  }
  return(param_list)
}


optim_wrapper <- function(dat, n_inits, param_ranges, model,
                          parallel = TRUE, n_cores = 3,
                          method = 'L-BFGS-B', maxit = 200,
                          trace = 0){
  
  # create param_list:
  if (model != 1 & model != 2 & model != 3) stop("Don't know which model to fit!")
  
  if (model == 1){
    init_list <- list(rmax = runif(n_inits, min = 1, max = 30),
                      p = runif(n_inits, min = param_ranges$p[1], max = param_ranges$p[2]),
                      q = runif(n_inits, min = param_ranges$q[1], max = param_ranges$q[2]),
                      z = matrix(runif(n_inits*length(levels(dat$sf_factor)), 
                                       min = param_ranges$z[1], max = param_ranges$z[2]), 
                                 ncol=length(levels(dat$sf_factor))))
    
    # need to mould init_list into a big matrix with params as columns:
    inits <- param_to_mat(init_list, optim_inits = TRUE)
    
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
    init_list <- list(rmax = runif(n_inits, min = 1, max = 30),
                      p = runif(n_inits, min = param_ranges$p[1], max = param_ranges$p[2]),
                      q = runif(n_inits, min = param_ranges$q[1], max = param_ranges$q[2]),
                      z = matrix(runif(n_inits*length(levels(dat$sf_factor)), 
                                       min = param_ranges$z[1], max = param_ranges$z[2]), 
                                 ncol=length(levels(dat$sf_factor))),
                      w = matrix(runif(n_inits*length(levels(dat$sf_factor)), 
                                       min = param_ranges$w[1], max = param_ranges$w[2]), 
                                 ncol=length(levels(dat$sf_factor))))
    
    # need to mould init_list into a big matrix with params as columns:
    inits <- param_to_mat(init_list, optim_inits = TRUE)
    
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
  
  
  
  fun <- function(i) {
    tryCatch({
      
      if (method == 'BFGS'){
        opt <- optim(par = inits[i,], log_lik, dat = dat, model = model,
                     method = 'BFGS',
                     control = list(maxit = maxit, fnscale = -1, trace = trace),
                     hessian = TRUE)        
      }
      
      if (method == 'L-BFGS-B'){
        opt <- optim(par = inits[i,], log_lik, dat = dat, model = model,
                     method = 'L-BFGS-B', lower = lower, upper = upper,
                     control = list(maxit = maxit, fnscale = -1, trace = trace),
                     hessian = TRUE)    
      }
      
      if (method == 'Nelder-Mead'){
        opt <- optim(par = inits[i,], log_lik, dat = dat, model = model,
                     control = list(maxit = maxit, fnscale = -1, trace = trace),
                     hessian = TRUE)  
      }
      
      return(opt)
      
    }, error = function(e){
      warning(conditionMessage(e), " when i was ", i)
      return(list(value = NA))
    })
  }
  
  if (parallel == TRUE){
    require(parallel)
    fit_out <- mclapply(1:NROW(inits),fun,
                        mc.preschedule = FALSE,
                        mc.cores = n_cores)  
  } else {
    fit_out <- lapply(1:NROW(inits),fun)  
  }
  
  # find best run out of all the initial conditions:
  vals <- NULL
  for (i in 1 : length(fit_out)){
    vals[i] <- fit_out[[i]]$value
  }
  
  ind <- which(vals == max(vals, na.rm = TRUE))
  
  if (length(ind) > 1) ind <- ind[1] # equal log likelihoods; pick first...
  if (length(ind) == 0) ind <- 1 # no iterations were able to produce useful estimates; take first NAs
  if (is.infinite(ind)) ind <- 1 # no iterations were able to produce useful estimates; take first NAs
  
  opt <- fit_out[[ind]]
  
  # also return the deviation in value between the different initialisations:
  opt$val_sd <- sd(vals, na.rm=TRUE)
  opt$n_NA <- length(vals[is.na(vals)==TRUE])
  
  # return best opt list:
  return(opt)  
}
