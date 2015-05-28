# fit simple nonlinear transducer model to data subset.

require(psybayes)
require(rstan)
require(plyr)

debug <- FALSE
iter <- 20000
n_saved_samples <- 5000

source(paste0(getwd(),"/funs/helpers.R"))

model_file <- paste0(getwd(),'/funs/single_level_glm.stan')
output_file <- paste0(getwd(),'/output/single_level_glm.RData')

load(paste0(getwd(),'/output/csf_data_subset.RData'))


# subset data for faster debugging of fits:
if (debug == TRUE){
  rand_rows <- sample(1:nrow(dat), size = 3000)
  dat <- dat[rand_rows,]
  iter <- 100
  n_saved_samples <- 100
}

# prepare data for stan models --------------------

X <- make_design_matrix(dat)

stan_dat <- list(N = nrow(dat), 
                 y = dat$correct,
                 X = X,
                 D = ncol(X),
                 ss = as.numeric(dat$subject),
                 S = length(levels(dat$subject))) 

# Load and fit stan model -------------------------
fit <- stan_sample(file_path=model_file, stan_dat, iter = iter, n_saved_samples = n_saved_samples)

print(fit)

save(fit, file=output_file)
