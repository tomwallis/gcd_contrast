## test scripts for helper functions.

rm(list = ls())

library(rstan)
library(ggplot2)
library(reshape2)

source(paste0(getwd(),'/funs/helpers.R'))

n_subjs <- 5
n_samples <- 50

# some test data (for plotting whole range):
test_dat <- expand.grid(subject = 1:n_subjs, alpha = seq(1.5, 6.5, length = 20), 
                        c_centre_target = seq(0, 0.25, length = 5))

test_dat$increment <- test_dat$c_centre_target * test_dat$alpha - test_dat$c_centre_target

####################################################################
# test transducer functions -------------------------------------
####################################################################

# make up param samples:
p <- matrix(rlnorm(n_samples*n_subjs, meanlog = log(2), sdlog = 2/30), ncol = n_subjs)
q <- matrix(rlnorm(n_samples*n_subjs, meanlog = log(0.4), sdlog = 0.4/30), ncol = n_subjs)
z <- matrix(rlnorm(n_samples*n_subjs, meanlog = log(0.1), sdlog = 0.1/30), ncol = n_subjs)
rmax <- matrix(rlnorm(n_samples*n_subjs, meanlog = log(6), sdlog = 10/30), ncol = n_subjs)

# perturb the individual subject samples a bit:
p[,1] <- p[,1] + 2
z[,2] <- z[,2] + 0.1
z[,3] <- exp(log(z[,3]) - 1)
q[,4] <- 0

# arrange as a list:
transducer_params <- list(p = p, q = q, z = z, rmax = rmax)

# generate predictions for test data:
test_dat$transducer_p <- trial_p_transducer(dat = test_dat, param_list = transducer_params)


# plot transducer functions -------------------------------------
fig <- ggplot(data = test_dat, aes(x = alpha, y = transducer_p)) +
  facet_grid(c_centre_target ~ subject) +
  geom_line() +
  ggtitle("P_hat for transducer")
fig


# contrast response function -------------------------------------
crf_dat <- expand.grid(subject = 1:n_subjs, contrast = exp(seq(-4, 0, length = 200)))
post_samples <- r_fun(param_list = transducer_params, c_vec = crf_dat$contrast, 
                      subj_ind = crf_dat$subject, priors = FALSE)

# post samples has one column per sample from the posterior. Plot each one:
post_samples <- melt(post_samples)
names(post_samples) <- c("trial", "sample", "R")
post_samples$subject <- crf_dat$subject
post_samples$contrast <- crf_dat$contrast


fig <- ggplot(data = post_samples, aes(x = contrast, y = R)) +
  facet_wrap(~ subject) +
  scale_x_log10() + 
  ggtitle("Contrast response functions")
for(samp in 1 : max(unique(post_samples$sample))){
  fig <- fig + geom_line(data = subset(post_samples, sample == samp), alpha = 0.2)
}
fig

# tvc functions -------------------------------------
tvc_samples <- tvc_fun(param_list = transducer_params, c_vec = crf_dat$contrast, 
                       subj_ind = crf_dat$subject, priors = FALSE)

# same as post samples, above:
tvc_samples <- melt(tvc_samples)
names(tvc_samples) <- c("trial", "sample", "Thresh")
tvc_samples$subject <- crf_dat$subject
tvc_samples$contrast <- crf_dat$contrast


fig <- ggplot(data = tvc_samples, aes(x = contrast, y = Thresh)) +
  facet_wrap(~ subject) +
  scale_x_log10() + scale_y_log10() + 
  ggtitle("TvC functions")
for(samp in 1 : max(unique(tvc_samples$sample))){
  fig <- fig + geom_line(data = subset(tvc_samples, sample == samp), alpha = 0.2)
}
fig

####################################################################
# test glm functions -------------------------------------
####################################################################

# make up param samples:
b0 <- matrix(rnorm(n_samples*n_subjs, mean = 0, sd = 0.1), ncol = n_subjs)
b1 <- matrix(rnorm(n_samples*n_subjs, mean = -0.5, sd = 0.1), ncol = n_subjs)
b2 <- matrix(rnorm(n_samples*n_subjs, mean = 1, sd = 0.2), ncol = n_subjs)
b3 <- matrix(rnorm(n_samples*n_subjs, mean = -0.2, sd = 0.01), ncol = n_subjs)

# perturb the individual subject samples a bit:
b1[,1] <- b1[,1] + 1
b2[,2] <- b2[,2] + 1
b3[,3] <- b3[,3] - 1
b3[,4] <- b3[,4] + 1

# arrange as array, then list:
beta <- array(c(b0, b1, b2, b3), dim = c(n_samples, n_subjs, 4))
glm_params <- list(beta = beta)

# generate predictions for test data:
test_dat$glm_p <- trial_p_glm(dat = test_dat, param_list = glm_params)

# plot glm functions -------------------------------------
fig <- ggplot(data = test_dat, aes(x = alpha, y = glm_p)) +
  facet_grid(c_centre_target ~ subject) +
  geom_line() +
  ggtitle("P_hat for GLM")
fig

# plot glm surface example -------------------------------------
surf_dat <- expand.grid(subject = 1:5,
                        c_centre_target = seq(0, 0.25, length = 100),
                        increment = seq(0, 0.75, length = 100))

# generate predictions:
surf_dat$glm_p <- trial_p_glm(dat = surf_dat, param_list = glm_params)

# plot:
fig <- ggplot(data = surf_dat, aes(x = c_centre_target, y = increment, z = glm_p)) +
  facet_wrap(~ subject) +
  geom_tile(aes(fill = glm_p)) +
  stat_contour() +
  ggtitle("GLM model surfaces")
fig

# wireframe plot:
library(lattice)
wireframe(glm_p ~ c_centre_target * increment, data = surf_dat)

