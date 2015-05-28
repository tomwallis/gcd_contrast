## fitting stan model on the NIP cluster...

iter <- 4000
n_saved_samples <- 4000

load(paste0(getwd(),'/output/simulated_data_3.RData'))
model_file <- paste0(getwd(),'/funs/gain_control_full.stan')
output_file <- paste0(getwd(),'/output/stan_sim_fits_3.RData')

norm_pool <- subset(dat, select = c_surround_band_0_target : c_surround_band_5_target)

stan_dat <- list(N = nrow(dat), 
                 y = dat$sim_correct,
                 ped = dat$c_centre_target,
                 ped_plus_inc = dat$c_centre_target + dat$increment,
                 ss = as.numeric(dat$subject),
                 S = length(levels(dat$subject)),
                 sf = as.numeric(dat$sf_factor),
                 SF = length(levels(dat$sf_factor)),
                 norm_pool = as.matrix(norm_pool)) 

# attach(stan_dat)
# stan_rdump(c('N', 'y', 'ped','ped_plus_inc','ss','S','sf','SF','norm_pool'), file = paste0(getwd(),'/output/simulated_data.txt'))
# detach(stan_dat)

# Load and fit stan model -------------------------
library(rstan)
library(psybayes)
fit <- stan_sample(file_path=model_file, stan_dat, iter = iter, n_saved_samples=n_saved_samples)

print(fit)

save(fit, file=output_file)  