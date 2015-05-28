## function and scripts to play with parameters of contrast response function.
# TSAW.

rm(list=ls())

library(ggplot2)
library(wutils)

# c is linear contrast (could be filter response)
# n is psychometric function slope
# z is the contrast corresponding to half the asymptotic performance
# rmax is the maximum response.

crf <- function(c,param_list){
  p <- param_list$p
  q <- param_list$q
  z <- param_list$z
  r <- (c^(p+q) / (c^p + z^p))
}    

crf_gc <- function(c,param_list){
  p <- param_list$p
  q <- param_list$q
  z <- param_list$z
  w <- param_list$w
  
  # denom <- TODO
  
  r <- (c^(p+q) / (c^p + z^p))
}    

source(paste(getwd(),'/funs/calc_tvc.R',sep=""))


# set up params:
p <- c(.5,2,4)
q <- c(.4, .2)
z <- c(.005,.02,.2, 1)
alpha_ratio <- c(1) # dprime increment to test for tvc.

# simulate some curves fit from this function to examine the effects of various params.
contrasts <- c(0, exp(seq(log(0.001),log(0.2),length=200))) # e.g. in percent.
gen <- expand.grid(c = contrasts, p = p, q = q, z = z, alpha_ratio = alpha_ratio)
param_list <- list(p = gen$p, q = gen$q, z = gen$z)
gen$r <- crf(c = gen$c, param_list = param_list)

fig <- ggplot(gen,aes(x=c,y=r,colour=factor(z))) + facet_grid(p ~ q) + geom_line()
fig <- fig + scale_x_log10(name="Contrast")
fig <- fig + ylab("Response (d' units)")
fig <- fig + ggtitle('Contrast response function (CRF)')
fig <- fig + theme(legend.position="top")
fig
ggsave(paste0(getwd(),'/demos/crf_params_demo.pdf'),width=6,height=5)

#--------------------
# TvC functions:

gen$thresh <- calc_tvc(crf_fun=crf,c=gen$c, param_list=param_list, d_prime_inc=gen$alpha_ratio[1])

fig <- ggplot(gen,aes(x=c,y=thresh,colour=factor(z))) + facet_grid(p ~ q) + geom_line()
fig <- fig + scale_x_log10(name="Contrast (%)")
fig <- fig + scale_y_log10("Threshold (% contrast)")
fig <- fig + ggtitle('Threshold versus Contrast (TvC)')
fig <- fig + ggtitle('Threshold versus contrast')
fig <- fig + theme(legend.position="top")
fig
ggsave(paste(getwd(),'/demos/tvc_params_demo.pdf',sep=""),width=6,height=5)




