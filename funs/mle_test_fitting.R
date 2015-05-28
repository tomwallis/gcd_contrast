## file to fit MLE for main data. Trying to keep script pretty bare
# to allow calling from cluster.

load(paste0(getwd(),'/output/csf_data_reduced.RData'))

source(paste0(getwd(),'/funs/mle_funs.R'))

n_inits <- 2
n_cores <- 2
maxit <- 10
method <- 'L-BFGS-B'

param_ranges <- list(p = c(0.2, 6),
                     q = c(0, 2),
                     z = c(0.0001, 20),
                     w = c(0, 1))

# fit each subject separately, seed inits in parallel ---------------

fit_frame <- data.frame()

for (subj in 1 : length(levels(dat$subject))){
  print(this_subject <- levels(dat$subject)[subj])
  sub_dat <- subset(dat, subject == this_subject)
  # remove extra levels of subject factor:
  sub_dat$subject <- factor(sub_dat$subject)
  
  # fit this subject's parameters:
  out <- optim_wrapper(dat = sub_dat, n_inits = n_inits, 
                       param_ranges = param_ranges, n_cores = n_cores,
                       method = method, maxit = maxit)
  
  # process fitted params to list:
  fitted_pars <- mat_to_param(out$par)
  
  this_frame <- expand.grid(subject = this_subject,
                            convergence = out$convergence,
                            value = out$value,
                            val_sd = out$val_sd,
                            n_NA = out$n_NA,
                            p = fitted_pars$p,
                            q = fitted_pars$q,
                            sf = sort(unique(dat$sf)))    
  
  this_frame$sf_factor <- levels(dat$sf_factor)
  this_frame$z <- fitted_pars$z
  this_frame$w <- fitted_pars$w
  
  fit_frame <- rbind(fit_frame, this_frame)
}


save(fit_frame, file = paste0(getwd(),'/output/mle_out_testing.RData'))
