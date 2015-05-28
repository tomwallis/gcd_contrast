## helper functions for plotting

# function to take a parameter set and output a data frame of samples for each param ----------
param_distribution <- function(fit, dat, model_num){
  require(rstan)
  params <- extract(fit)
  n_samples <- dim(params[[1]])[1]
  
  types <- c('Prior','Posterior')
  
  if(model_num == 1){
    model <- "GLM"
    coeffs <- c('beta_1','beta_2','beta_3')
  }
  
  if(model_num == 2){
    model <- "Multi GLM"
    coeffs <- c('beta_1','beta_2','beta_3')
  }
  
  if(model_num == 3){
    model <- "Transducer A"
    coeffs <- c('p','q','z','rmax')
  }
  
  if(model_num == 4){
    model <- "Transducer B"
    coeffs <- c('p','q','z','rmax')
  }
  
  model_samps <- expand.grid(subject = unique(dat$subject), model = model, type = types,
                             coeff = coeffs, sample = 1:n_samples)
  model_samps$value <- NA
  
  for (type in 1:length(types)){
    this_type <- types[type]
    for (coeff in 1:length(coeffs)){
      this_coeff <- coeffs[coeff]
      for (subj in 1 : length(levels(model_samps$subject))){
        this_subj <- levels(model_samps$subject)[subj]
        
        this_row <- model_samps$subject == this_subj & 
          model_samps$type == this_type &
          model_samps$coeff == this_coeff
        
        if(model_num == 1 | model_num == 2){
          if (this_type == 'Prior'){
            this_par <- "prior_beta"
          } else {
            this_par <- "beta"
          }
          text = paste0("this_vec <- params$",this_par,"[,", subj, ", ", coeff, "]")
          eval(parse(text = text))
        }
        
        if(model_num == 3 | model_num == 4){
          if (this_type == 'Prior'){
            this_par <- paste0("prior_",this_coeff)
          } else {
            this_par <- this_coeff
          }
          text = paste0("this_vec <- params$",this_par,"[,", subj, "]")
          eval(parse(text = text))
        }
        
        # place into model samps:
        model_samps$value[this_row] <- this_vec
      }
    }
  }
  return(model_samps)
}


# generate model predictions for a surface ----------

surface_prediction <- function(fit, dat, model_num, subject = 1){
  require(rstan)
  params <- extract(fit)
  n_samples <- dim(params[[1]])[1]
  
  # values for pedestal and increment picked based on data range:
#   surf_dat <- expand.grid(subject = subject,
#                           c_centre_target = seq(0, 0.1, length = 100),
#                           increment = seq(0, 0.5, length = 100))
  surf_dat <- expand.grid(subject = subject,
                        c_centre_target = exp(seq(-9, -0.5, length = 25)),
                        increment = exp(seq(-9, -0.5, length = 25)))
  
  # generate predictions:
  if(model_num == 1){
    surf_dat$model <- "GLM"
    surf_dat$p <- trial_p_glm(dat = surf_dat, param_list = params)
  }
  
  if(model_num == 2){
    surf_dat$model <- "Multi GLM"
    surf_dat$p <- trial_p_glm(dat = surf_dat, param_list = params)
  }
  
  if(model_num == 3){
    surf_dat$model <- "Transducer A"
    surf_dat$p <- trial_p_transducer(dat = surf_dat, param_list = params)
  }
  
  if(model_num == 4){
    surf_dat$model <- "Transducer B"
    surf_dat$p <- trial_p_transducer(dat = surf_dat, param_list = params)
  }
  
  return(surf_dat)
}

# generate model predictions for a surface ----------

alpha_prediction <- function(fit, dat, model_num, breaks){
  require(rstan)
  require(dplyr)
  params <- extract(fit)
  n_samples <- dim(params[[1]])[1]
  
  # for each trial of the data frame, return predictions:
  alpha_dat <- subset(dat, select = c("subject", "alpha", "c_centre_target", "increment"))
  
  # generate predictions:
  if(model_num == 1){
    alpha_dat$model <- "GLM"
    phat_matrix <- phat_samples_glm(dat = alpha_dat, param_list = params)
  }
  
  if(model_num == 2){
    alpha_dat$model <- "Multi GLM"
    phat_matrix <- phat_samples_glm(dat = alpha_dat, param_list = params)
  }
  
  if(model_num == 3){
    alpha_dat$model <- "Transducer A"
    phat_matrix<- phat_samples_transducer(dat = alpha_dat, param_list = params)
  }
  
  if(model_num == 4){
    alpha_dat$model <- "Transducer B"
    phat_matrix <- phat_samples_transducer(dat = alpha_dat, param_list = params)
  }
  
  alpha_dat$c_centre_target_factor <- cut(alpha_dat$c_centre_target,breaks=equalBreaks,include.lowest=TRUE)
  levels(alpha_dat$c_centre_target_factor) <- c('Pedestal low','Pedestal high')
  
  # for each subject, alpha level and pedestal factor, calculate hdi on predictions:
  summary_dat <- expand.grid(subject = levels(alpha_dat$subject),
                             alpha = unique(alpha_dat$alpha),
                             c_centre_target_factor = levels(alpha_dat$c_centre_target_factor),
                             model = factor(unique(alpha_dat$model)))
  summary_dat$y <- NA
  summary_dat$ymin <- NA
  summary_dat$ymax <- NA

  for (subj in 1 : length(levels(summary_dat$subject))){
    for (alpha in 1 : length(unique(summary_dat$alpha))){
      for (ped in 1 : length(levels(summary_dat$c_centre_target_factor))){
        this_subj <- levels(summary_dat$subject)[subj]
        this_alpha <- unique(summary_dat$alpha)[alpha]
        this_ped <- levels(summary_dat$c_centre_target_factor)[ped]
        this_rows_alpha <- which(alpha_dat$subject == this_subj &
                             alpha_dat$alpha == this_alpha &
                             alpha_dat$c_centre_target_factor == this_ped)

        this_rows <- which(summary_dat$subject == this_subj &
                             summary_dat$alpha == this_alpha &
                             summary_dat$c_centre_target_factor == this_ped)
        
        # all the p_hats for this alpha, subject and pedestal:
        this_samples <- as.vector(phat_matrix[this_rows_alpha, ])

        if (length(this_samples > 0)){
          summary_dat$y[this_rows] <- mean(this_samples)
          hdis <- hdi(this_samples)
          summary_dat$ymin[this_rows] <- hdis[1]
          summary_dat$ymax[this_rows] <- hdis[2]          
        }
      }
    }
  }
  
  return(summary_dat)
}


