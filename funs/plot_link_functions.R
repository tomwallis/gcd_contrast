# compare link functions (d prime and my weibull approx):

library(ggplot2)

# the d' range to test:
d_range <- seq(0,5,l=200)

# Normal SDT --------------------------

# from the psyphy function "dprime.mAFC":
pc_sdt <- function(dprime,m){
  pr <- function(x) dnorm(x - dprime) * pnorm(x)^(m - 1)
  return(integrate(pr, lower = -Inf, upper = Inf)$value)
}

ah_approx <- function(dprime, m){
  pnorm(dprime / sqrt(2))^log2(m)
}

wb_fun <- function(dprime, m, shape = 1, scale = 1){
  raw_pc <- pweibull(dprime, shape = shape, scale = scale)
  gamma <- 1 / m
  return ( gamma + (1 - gamma) * raw_pc )
}

join <- function(i, dprime, m, name, shape = 1, scale = 1){
  switch(name[i], 
         SDT = pc_sdt(dprime[i], m[i]),
         ah = ah_approx(dprime[i], m[i]),
         Weibull = wb_fun(dprime[i], m[i], shape = shape[i], scale = scale[i]))
}


# Generate predictions ----------------

dat <- expand.grid(dprime = d_range, m = 4, shape = 1.481270, scale = 1.545903,
                   name = c('SDT','Weibull'))

dat$y <- sapply(1:nrow(dat),join, dprime = dat$dprime, m = dat$m, name = dat$name, shape = dat$shape, scale = dat$scale)


library(ggplot2)
fig <- ggplot(dat, aes(x = dprime, y = y, colour = name)) + geom_line() + 
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'top') + 
  ylim(0,1) + 
  ylab('p(correct)') + xlab('Delta R') + 
  scale_colour_brewer(name = "", type = "qual", palette = 2)

ggsave(fig, file=paste0(getwd(),'/figs/link_functions.pdf'), width = 3, height = 3)

