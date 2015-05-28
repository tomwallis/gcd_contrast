## fitting stan model on the NIP cluster...

require(psybayes)
require(rstan)
require(plyr)

debug <- FALSE
iter <- 30000
n_saved_samples <- 10000

model_file <- paste0(getwd(),'/funs/gain_control_full.stan')
output_file <- paste0(getwd(),'/output/gain_control_full_long.RData')

load(paste0(getwd(),'/output/csf_data_reduced.RData'))


# subset data for faster debugging of fits:
if (debug == TRUE){
  rand_rows <- sample(1:nrow(dat), size = 3000)
  dat <- dat[rand_rows,]    
}

# prepare data for stan models --------------------
norm_pool <- subset(dat, select = c_surround_band_0_target : c_surround_band_5_target)

stan_dat <- list(N = nrow(dat), 
                 y = dat$correct,
                 ped = dat$c_centre_target,
                 ped_plus_inc = dat$c_centre_target + dat$increment,
                 ss = as.numeric(dat$subject),
                 S = length(levels(dat$subject)),
                 sf = as.numeric(dat$sf_factor),
                 SF = length(levels(dat$sf_factor)),
                 norm_pool = as.matrix(norm_pool)) 

# Load and fit stan model -------------------------
fit <- stan_sample(file_path=model_file, stan_dat, iter = iter, n_saved_samples = n_saved_samples)

print(fit)

save(fit, file=output_file)
