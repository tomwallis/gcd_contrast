## file to fit MLE, seeding from generating values of simulated data. Normalised deviance should be close to 1.

rm(list=ls())

load(paste0(getwd(),'/output/simulated_data_1.RData'))

# change dat$correct to sim_correct (don't have to dig through mle_funs and put conditional flags):
dat$correct <- dat$sim_correct

source(paste0(getwd(),'/funs/mle_funs.R'))

param_ranges <- list(rmax = c(0.001, 30),
                     p = c(0.2, 6),
                     q = c(0, 0.6),
                     z = c(0.0001, 1),
                     w = c(0, 1))

subj <- 3
print(this_subject <- levels(dat$subject)[subj])
sub_dat <- subset(dat, subject == this_subject)
# remove extra levels of subject factor:
sub_dat$subject <- factor(sub_dat$subject)

# use initial values for subject 2 in simulated data:
inits <- param_to_mat(gen_params)
inits <- inits[subj,]

initial_val <- log_lik(par = inits, dat = sub_dat, model = 1)

# set up lower and upper bounds for params:
lower <- c(param_ranges$rmax[1],
           param_ranges$p[1], 
           param_ranges$q[1],
           rep(param_ranges$z[1], times=6))
upper <- c(param_ranges$rmax[2],
           param_ranges$p[2], 
           param_ranges$q[2],
           rep(param_ranges$z[2], times=6))


opt <- optim(par = inits, log_lik, dat = sub_dat, model = 1,
             method = 'L-BFGS-B', lower = lower, upper = upper,
             control = list(maxit = 100, trace = 6, fnscale = -1, factr = 1e-5))    

save(opt, initial_val, file = paste0(getwd(),'/output/mle_test_single_out.RData'))

