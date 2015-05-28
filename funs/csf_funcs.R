## Functions to produce CSFs given some parameters. Default parameters approximated from Watson & Ahumada, Table 5.
#
# TSAW

#----------------------
# Yang-Qi-Makous CSF (same number of parameters as log-parabola but continuous function):
yqm_csf <- function(f=exp(seq(log(.5),log(30),length=50)),
                    f_zero=7.06,
                    gain=466,
                    f_one=0.69,
                    a=7.77){
  
  # output sensitivity for a range of frequencies f.
  # params are 
  s <- exp(-f / f_zero) / (1 + (a / (1 + (f / f_one)^2) ))
  # scale by gain:
  s <- s * gain
  return(s)
}

#----------------------
# Log_parabola CSF (like is used by the qCSF):
lp_csf <- function(f=exp(seq(log(.5),log(30),length=50)),
                   f_max=3.23,
                   gain=214.46,
                   beta=2.49,
                   a=0.71){
  # log parabola CSF with low frequency truncation.
  
  # output sensitivity for a range of frequencies f.
  # params are freq (log), f_max (log units), gain (log units), beta (log units), truncation (log)
  tau_decay <- 2
  kappa <- log10(tau_decay)
  beta_prime <- log10(2 * beta)
  s_prime <- log10(gain) - kappa * ( (log10(f) - log10(f_max)) / (beta_prime / 2) )^2
  s <- s_prime
  s[f < f_max & s_prime < (log10(gain) - a)] <- log10(gain) - a
  lin_s <- 10^s
  return(lin_s)
}


#------------------
# function to convert sensitivity output by a csf into the z parameter in a transducer function (see Haun and Essock, equation B5):
s_to_z <- function(s,d_prime=1,p=2,q=0.4,rmax=6){
  # default values taken from population parameters in 4 parameter model.
  z <- (1/s) * ((rmax / d_prime)* (1/s)^q - 1)^(1/p)
  return(z)
}

s_to_z_vec <- function(i,s,d_prime=1,p=c(0.5,2),q=c(0.3,0.4),rmax=c(6,6)){
  # default values taken from population parameters in 4 parameter model.
  z <- (1/s[i]) * ((rmax[i] / d_prime)* (1/s[i])^q[i] - 1)^(1/p[i])
  return(z)
}

#------------------
# # plot examples. for more see explore scripts in development files.
# library(ggplot2)
# 
# contrast <- exp(seq(log(.5),log(30),length=100))
# dat <- expand.grid(f=contrast,model=c('lp','yqm'))
# dat$y <- NA
# dat$y[dat$model=='lp'] <- lp_csf(f=contrast)
# dat$y[dat$model=='yqm'] <- yqm_csf(f=contrast)
# 
# fig <- ggplot(dat,aes(x=f,y=y,colour=model))
# fig <- fig + geom_line()
# fig <- fig + scale_x_log10()
# fig

