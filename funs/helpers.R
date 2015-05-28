## generic helper functions to be sourced for running analysis.
# TSAW

log_convert <- function(x){
  x[x==0] <- min(x[x > 0]) # set to minimum nonzero value
  x <- log(x)
  return(x)
}

############################################################
# Transducer functions
############################################################

# return p correct for each trial (vectorised) -----------
trial_p_transducer <- function(dat, param_list, priors = FALSE, return_range = FALSE){
  phat_matrix <- phat_samples_transducer(dat = dat, 
                                  param_list = param_list, 
                                  priors = priors)
  
  # summarise the phat matrix somehow:
  if(return_range == TRUE){
    y <- rowMeans(phat_matrix)
    library(psybayes)
    hdis <- col_hdi(t(phat_matrix))
    ymin <- hdis[1,]
    ymax <- hdis[2,]
    
    return(list(y = y, ymin = ymin, ymax = ymax))
    
  } else {
    return(rowMeans(phat_matrix))
  }
}


phat_samples_transducer <- function(dat, param_list, priors = FALSE){
  ped <- dat$c_centre_target
  ped_plus_inc <- dat$c_centre_target + dat$increment
  subj <- as.numeric(dat$subject)
  
  delta_r <- delta_r_fun(param_list = param_list, 
                         ped = ped, ped_plus_inc = ped_plus_inc, 
                         subj_ind = subj,
                         priors = priors)
  
  phat <- link_fun(delta_r)
  # p is a matrix with rows as trials from dat and cols as samples from mcmc.
  return(phat)
}

# compute delta r --------------
delta_r_fun <- function(param_list, ped, ped_plus_inc, subj_ind, priors){
  # compute response to ped:
  r_ped <- r_fun(param_list = param_list, c_vec = ped, 
                 subj_ind = subj_ind,
                 priors = priors)
  
  # compute response to ped + inc:
  r_ped_plus_inc <- r_fun(param_list = param_list, c_vec = ped_plus_inc, 
                          subj_ind = subj_ind, 
                          priors = priors)
  
  # difference:
  delta_r <- r_ped_plus_inc - r_ped
  
  # check for zero:
  delta_r <- pmax(delta_r, 0) # pmax should work for scalar, vector or matrix.
  return(delta_r)
}


# contrast reponse function ------------------
r_fun <- function(param_list, c_vec, subj_ind, priors){
  
  if(priors == TRUE){
    p <- param_list$prior_p[, subj_ind]
    q <- param_list$prior_q[, subj_ind]
    z <- param_list$prior_z[, subj_ind]
    rmax <- param_list$prior_rmax[, subj_ind]
  } else {
    p <- param_list$p[, subj_ind]
    q <- param_list$q[, subj_ind]
    z <- param_list$z[, subj_ind]
    rmax <- param_list$rmax[, subj_ind]    
  }
  
  # transpose to have rows == trials and cols == samples:
  p <- t(p)
  q <- t(q)
  z <- t(z)
  rmax <- t(rmax)
  
  denom <- z^p + c_vec^p
  numer <- rmax * c_vec^(p + q)
  
  r <- numer / denom
  return(r)  
} 


# link function, returning p(correct) ------------------------------
link_fun <- function(delta_r, m = 4){
  raw_pc <- pweibull(delta_r, shape = 1.481270, scale = 1.545903)
  gamma <- 1 / m
  return ( gamma + (1 - gamma) * raw_pc )
}

# tvc approximation ------------------
tvc_fun <- function(param_list, c_vec=exp(seq(-4,0,length=100)), subj_ind, priors = FALSE,
                     d_prime_inc=1, grad_step=.001){
  
  r_standard <- r_fun(param_list, c_vec, subj_ind = subj_ind, priors = priors)
  r_increment <- r_fun(param_list, (c_vec + grad_step), subj_ind = subj_ind, priors = priors)
  
  delta_r <- r_increment - r_standard
  
  # threshold = derivative of contrast with respect to response:
  thresh <- d_prime_inc * (grad_step / delta_r)
  return(thresh)
}

############################################################
# GLM helper functions
############################################################

make_design_matrix <- function(dat, return_scale_params = FALSE, scale_factors = NULL){
  # fit to log contrast and increment: first set zeros to lowest non-zero val:
  dat$c_centre_target[dat$c_centre_target==0] <- min(dat$c_centre_target[dat$c_centre_target>0])
  dat$increment[dat$increment==0] <- min(dat$increment[dat$increment>0])
  
  # design matrix:
  X <- model.matrix(~ log(dat$c_centre_target) + log(dat$increment))
  
  # normalise:
  if (is.null(scale_factors)){
    scaled <- scale(X[, -1])
  } else {
    scaled <- scale(X[, -1], center = scale_factors$center, scale = scale_factors$scale)
  }
  
  X[, -1] <- scaled
  
  if (return_scale_params == TRUE){
    return(list(X = X, center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale")))
  } else {
    return(X)
  }
  
}

trial_p_glm <- function(dat, param_list, priors = FALSE, 
                        return_range = FALSE, full_model = FALSE,
                        scale_factors = NULL){  
  phat_matrix <- phat_samples_glm(dat = dat, 
                                  param_list = param_list, 
                                  priors = priors,
                                  full_model = full_model,
                                  scale_factors = scale_factors)
  
  # summarise the phat matrix somehow:
  if(return_range == TRUE){
    y <- rowMeans(phat_matrix)
    library(psybayes)
    hdis <- col_hdi(t(phat_matrix))
    ymin <- hdis[1,]
    ymax <- hdis[2,]
    
    return(list(y = y, ymin = ymin, ymax = ymax))
    
  } else {
    return(rowMeans(phat_matrix))
  }
}

phat_samples_glm <- function(dat, param_list, priors = FALSE,
                             full_model = FALSE,
                             scale_factors = NULL){  
  pred_p <- function(beta_vec){
    eta <- rowSums(X * beta_vec)
    p <- 0.25 + (1 - 0.25) * plogis(eta)
    return(p)
  }
  
  by_sample <- function(i){
    this_beta <- betas[i, , ]
    this_p <- pred_p(this_beta)
  }
  
  if (full_model == TRUE){
    if (is.null(scale_factors)){
      X <- make_design_matrix_full(dat)
    } else {
      X <- make_design_matrix_full(dat, scale_factors = scale_factors)
    }
    
  } else {
    if (is.null(scale_factors)){
      X <- make_design_matrix(dat)  
    } else {
      X <- make_design_matrix(dat, scale_factors = scale_factors)  
    }
  }

  subj_ind <- as.numeric(dat$subject)
  n_samples <- dim(param_list$beta)[1]
  
  # find beta matrix for all subjects:
  if (priors == FALSE){
    betas <- param_list$beta[, subj_ind, ]
  } else {
    betas <- param_list$prior_beta[, subj_ind, ]
  }
  
  p_hat <- sapply(1:n_samples, by_sample)
  
  # phat is a matrix where rows are the same as dat, and columns are the predictions for each mcmc sample.
  return(p_hat)
}

############################################################
# Full GLM helper functions
############################################################

make_design_matrix_full <- function(dat, return_scale_params = FALSE, scale_factors = NULL){

  # calling this from surface plotting functions means you have log_K3D etc already in the frame. Test.
  if (is.null(dat$log_K3D)){
    # do log conversions:
    dat$log_K3D <- log_convert(dat$K_3D_21_target)
    dat$log_ped <- log_convert(dat$c_centre_target)
    dat$log_inc <- log_convert(dat$increment)
    dat$log_em_cum <- log_convert(dat$em_cumDist)
  }
  
  # design matrix:
  X <- model.matrix(~ log_ped + log_inc + log_K3D + 
                      log_em_cum + sf_factor, data = dat)

  scale_cols <- 2:5 # this is quite an ugly hack, and will need to change if the model does.
  # normalise:
  if (is.null(scale_factors)){
    scaled <- scale(X[, scale_cols])
  } else {
    scaled <- scale(X[, scale_cols], center = scale_factors$center, scale = scale_factors$scale)
  }
  
  X[, scale_cols] <- scaled
  
  if (return_scale_params == TRUE){
    return(list(X = X, center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale")))
  } else {
    return(X)
  }
  
}

