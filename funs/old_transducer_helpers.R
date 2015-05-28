## generic helper functions to be sourced for running analysis.
# TSAW

# link function, returning p(correct) ------------------------------

link_fun <- function(delta_r, m = 4){
  raw_pc <- pweibull(delta_r, shape = 1.481270, scale = 1.545903)
  gamma <- 1 / m
  return ( gamma + (1 - gamma) * raw_pc )
}

######################################################################
# model input / output functions -----------------------------------
######################################################################

# function to calculate p correct through each trial (vectorised) -----------
trial_p <- function(dat, param_list, norm_pool, model, return_delta_r = FALSE){
  ped <- dat$c_centre_target
  ped_plus_inc <- dat$c_centre_target + dat$increment
  subj <- as.numeric(dat$subject)
  targ_sf <- as.numeric(dat$sf_factor)
  
  
  # check if you've passed a subject's data as a subset:
  #   if(length(unique(dat$subject))==1) subj <- rep(1, l=NROW(ped))
  
  delta_r <- delta_r_fun(param_list = param_list, 
                         ped = ped, ped_plus_inc = ped_plus_inc, 
                         norm_pool = norm_pool, 
                         targ_sf = targ_sf, subj_ind = subj,
                         model = model)
  
  p <- link_fun(delta_r)
  
  # If p is a matrix returned by stan (i.e. trials by samples, summarise)
  if (is.matrix(p) & NCOL(p)>1) browser()
  
  if (return_delta_r == TRUE){
    return(list(p = p, delta_r = delta_r))    
  } else {
    return(p)  
  }
}

# compute delta r --------------
delta_r_fun <- function(param_list, ped, ped_plus_inc, norm_pool, targ_sf, subj_ind, model){
  # compute response to ped:
  r_ped <- r_fun(param_list = param_list, contrast = ped, 
                 norm_pool = norm_pool, 
                 targ_sf = targ_sf, subj_ind = subj_ind,
                 model = model)
  
  # compute response to ped + inc:
  r_ped_plus_inc <- r_fun(param_list = param_list, contrast = ped_plus_inc, 
                          norm_pool = norm_pool, 
                          targ_sf = targ_sf, subj_ind = subj_ind,
                          model = model)
  
  # difference:
  delta_r <- r_ped_plus_inc - r_ped
  
  # check for zero:
  delta_r <- pmax(delta_r, 0) # pmax should work for scalar, vector or matrix.
  
  return(delta_r)
}


# contrast reponse function ------------------
# needs to be happy taking input from single rows, as columns, or as stan samples (3D array).
r_fun <- function(param_list, contrast, norm_pool, targ_sf, subj_ind, model){
  if (model != 1 & model != 2 & model != 3) stop("Don't know which model to fit!")
  
  # check input dimensions to decide what we're trying to do:
  if (is.data.frame(param_list)==TRUE){
    # param_list is a labelled data frame. Probably called from generating value.

    if (length(param_list$p) != 1) browser() # something is wrong!
    
    numerator <- param_list$rmax * contrast^(param_list$p + param_list$q)
    
    ### model-dependent code ###
    if (model == 1){
      # nonlinear transducer model.
      text <- paste0('this_z <- param_list$z',targ_sf)
      eval(parse(text=text))
      denominator <- this_z^param_list$p + contrast^param_list$p
    }
    
    if (model == 2 | model == 3){
      # gain control model (vectorised).
      text <- paste0('this_z <- param_list$z',targ_sf)
      eval(parse(text=text))
      this_w <- as.numeric(subset(param_list, select = w1:w6))
      w_norms <- norm_pool^param_list$p %*% this_w
      denominator <- this_z^param_list$p + w_norms
    }
    ### /model-dependent code ###
    return(r <- numerator / denominator )
  } 
  
  if(is.list(param_list) & length(param_list$p) == 1 & length(param_list$z) == 6){
    # calling from optim, for one sample of one subject.
    numerator <- param_list$rmax * contrast^(param_list$p + param_list$q)
    
    ### model-dependent code ###
    if (model == 1){
      # nonlinear transducer model.
      this_z <- as.numeric(param_list$z[targ_sf])
      denominator <- this_z^param_list$p + contrast^param_list$p
    }
    
    if (model == 2 | model == 3){
      # gain control model (vectorised).
      this_z <- as.numeric(param_list$z[targ_sf])
      this_w <- as.numeric(param_list$w)
      w_norms <- rowSums(norm_pool^param_list$p * this_w)
      denominator <- this_z^param_list$p + w_norms
      denominator <- as.vector(denominator)
    }
    ### /model-dependent code ###
    return(r <- numerator / denominator )
  }
  
  if(length(dim(param_list$z))==2 & is.list(param_list)){
    # calling from generating function -- z is a matrix with rows = sf and cols = subj.
    
    #     if (length(param_list$p) != 1) browser() # something is wrong!
    
    numerator <- param_list$rmax[subj_ind] * 
      contrast^(param_list$p[subj_ind] + param_list$q[subj_ind])
    
    
    # need to index z properly... can combine over models.
    z <- param_list$z[ cbind(targ_sf, subj_ind) ]
    
    
    ### model-dependent code ###
    if (model == 1){
      # nonlinear transducer model.
      denominator <- z^param_list$p[subj_ind] + contrast^param_list$p[subj_ind]
    }
    
    if (model == 2 | model == 3){
      # gain control model (vectorised).
      w_mat <- t(param_list$w[, subj_ind])
      w_norms <- rowSums(norm_pool^param_list$p[subj_ind] * w_mat)
      denominator <- z^param_list$p[subj_ind] + w_norms
    }
    return(r <- numerator / denominator )
  }
  
  if(length(dim(param_list$z))==3 & is.list(param_list)){
    # processing stan fit -- z is a 3D array with 
    # dim[1] = samples,
    # dim[2] = sf
    # dim[3] = subj.
    
    n_samples <- dim(param_list$p)[1]
    
    samp_fun <- function(i){
      # for each mcmc sample...
      numerator <- param_list$rmax[i, subj_ind] * 
        contrast^(param_list$p[i, subj_ind] + param_list$q[i, subj_ind])
      
      # need to index z properly... can combine over models.
      z <- param_list$z[i, cbind(targ_sf, subj_ind) ]
      
      ### model-dependent code ###
      if (model == 1){
        # nonlinear transducer model.
        denominator <- z^param_list$p[i, subj_ind] + contrast^param_list$p[i, subj_ind]
      }
      
      if (model == 2 | model == 3){
        # gain control model (vectorised).
        w_mat <- t(param_list$w[i, , subj_ind])
        w_norms <- rowSums(norm_pool^param_list$p[i, subj_ind] * w_mat)
        browser() # check that this is doing the right thing.
        denominator <- z^param_list$p[i, subj_ind] + w_norms
      }
      return(r <- numerator / denominator )
    }
    r <- sapply(1:n_samples, FUN = samp_fun)
    browser() # check dimensions of r, how to pass it to other functions.
    # TODO: parallelise here?
    return(r)
  }
  
} 

######################################################################
######################################################################

# #-------------------
# # function to turn large data frame into x and y, for contrast response fitting and other variables of interest.
# process_df_to_xy <- function(data_frame, model_name='default', design = NA, na_omit = FALSE, split_factor=NA){
#   
#   # take sf_factor out of data frame:
#   if (is.na(split_factor)){
#     x_vars <- c('subject','alpha','c_centre_target','sf','sf_factor','increment')
#     x <- data_frame[,x_vars]
#     y <- data_frame[,'correct']    
#   } else {
#     x_vars <- c('subject','alpha','c_centre_target','sf','sf_factor','increment',split_factor)
#     x <- data_frame[,x_vars]
#     y <- data_frame[,'correct']          
#   }  
#   
#   # omit nas if applicable.
#   xy <- cbind(x,y)
#   if(na_omit==TRUE) xy <- na.omit(xy)
#   
#   x <- subset(xy,select = -y)
#   y <- xy$y
#   
#   # apply model matrix:
#   if(any(is.na(design))==FALSE) {
#     x <- model.matrix(design,data=x)
#   }
#   
#   return(list(x=x,y=y))
# }

