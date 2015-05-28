## testing different "link" functions.
# different ways to convert the response difference
# to probability correct.

rm(list=ls())

library(wutils)
library(ggplot2)

# the d' range to test:
d_range <- seq(0,4,l=400)

# Normal SDT --------------------------

# from the psyphy function "dprime.mAFC":
pc_sdt <- function(dprime,m){
  pr <- function(x) dnorm(x - dprime) * pnorm(x)^(m - 1)
  return(integrate(pr, lower = -Inf, upper = Inf)$value)
}

# Andrew Haun analytic approximation -------------

ah_approx <- function(dprime, m){
  pnorm(dprime / sqrt(2))^log2(m)
}

# same as above with logit inverse instead of pnorm (faster for stan) -------------

ah_logit <- function(dprime, m){
  plogis(dprime / sqrt(2))^log2(m)
}

# Weibull link function, rescaled with lower asymptote -------------

wb_fun <- function(dprime, m, shape = 1, scale = 1){
  raw_pc <- pweibull(dprime, shape = shape, scale = scale)
  gamma <- 1 / m
  return ( gamma + (1 - gamma) * raw_pc )
}

# manually calculate weibull and check it's the same.
cum_weib <- function(x, shape = 1, scale = 1){
  1 - exp( - (x / scale)^shape)
}


# Joining function ---------------------

join <- function(i, dprime, m, name, shape = 1, scale = 1){
  switch(name[i], 
         sdt = pc_sdt(dprime[i], m[i]),
         ah = ah_approx(dprime[i], m[i]),
         logit = ah_logit(dprime[i], m[i]),
         weibull = wb_fun(dprime[i], m[i], shape = shape[i], scale = scale[i]))
}


# Generate predictions ----------------

dat <- expand.grid(dprime = d_range, m = c(2, 4, 8), shape = c(1, 2, 1.481270), scale = 1.545903,
                   name = c('sdt','ah','logit','weibull'))

dat$y <- sapply(1:nrow(dat),join, dprime = dat$dprime, m = dat$m, name = dat$name, shape = dat$shape, scale = dat$scale)

           
# plot ----------------
fig <- ggplot(dat, aes(x = dprime, y = y, colour = name)) + geom_line() + facet_grid(m ~ shape)
fig


# optimise the parameters of the weibull function to fit the SDT prediction ----------------

# function to difference values for all x:
diff_fun <- function(par, x = d_range, m = 4){
  
  vec_sdt <- function(i, x, m) pc_sdt(x[i], m)
  vec_weib <- function(i, x, m, shape, scale) wb_fun(x[i], m, shape = shape, scale = scale)
  
  y <- sapply(1:length(x), vec_sdt, x = x, m = m)
  y_hat <- sapply(1:length(x), vec_weib, x = x, m = m, shape = par[1], scale = par[2])

  # compute squared errors:
  er <- (y - y_hat)^2
  return(sum(er))  
}

par_vec <- function(n_inits){
  
  inits <- matrix(runif(n_inits*2, min = 0, max = 20), ncol = 2)
  
  pars <- matrix(rep(NA,l=2 * n_inits),ncol=2)
  val <- rep(NA, l = n_inits)
  
  for(i in 1 : n_inits){
    opt <- optim(par = inits[i,], diff_fun)
    pars[i,] <- opt$par
    val[i] <- opt$value
  }
  
  return(list(par = pars, value = val))
}


# opt_out <- par_vec(50)
# opt_ind <- which(opt_out$value == min(opt_out$value))
# print(opt_out$par[opt_ind,])

