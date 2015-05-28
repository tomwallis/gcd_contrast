## file to fit MLE for main data. Trying to keep script pretty bare
# to allow calling from cluster.

load(paste0(getwd(),'/output/csf_data_reduced.RData'))

library(hydroPSO)
source(paste0(getwd(),'/funs/mle_funs.R'))

n_cores <- 5
maxit <- 2000

param_ranges <- list(p = c(0.2, 6),
                     q = c(0, 2),
                     z = c(0.0001, 20),
                     w = c(0, 1))

# set up lower and upper bounds for params:
lower <- c(param_ranges$p[1], 
           param_ranges$q[1],
           rep(param_ranges$z[1], times=6),
           rep(param_ranges$w[1], times=6))
upper <- c(param_ranges$p[2], 
           param_ranges$q[2],
           rep(param_ranges$z[2], times=6),
           rep(param_ranges$w[2], times=6))

# fit each subject separately, seed inits in parallel ---------------
fit_frame <- data.frame()

for (subj in 1 : length(levels(dat$subject))){
  print(this_subject <- levels(dat$subject)[subj])
  sub_dat <- subset(dat, subject == this_subject)
  # remove extra levels of subject factor:
  sub_dat$subject <- factor(sub_dat$subject)
  
  # fit this subject's parameters:
  out <- hydroPSO(fn = log_lik, dat = sub_dat,
                  lower = lower, upper = upper,
                  control = list(maxit = 1000, MinMax = 'max', reltol = 1e-8, 
                                 write2disk=FALSE, normalise = TRUE,
                                 parallel = 'parallel', par.nnodes = n_cores))    
  
  # process fitted params to list:
  fitted_pars <- mat_to_param(out$par)
  
  this_frame <- expand.grid(subject = this_subject,
                            convergence = out$convergence,
                            message = out$message,
                            value = out$value,
                            p = fitted_pars$p,
                            q = fitted_pars$q,
                            sf = sort(unique(dat$sf)))    
  
  this_frame$sf_factor <- levels(dat$sf_factor)
  this_frame$z <- fitted_pars$z
  this_frame$w <- fitted_pars$w
  
  fit_frame <- rbind(fit_frame, this_frame)
}

save(fit_frame, file = paste0(getwd(),'/output/swarm_main_out.RData'))
