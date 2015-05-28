## function and scripts to play with parameters of contrast response function.
# TSAW.

rm(list=ls())

library(ggplot2)
library(wutils)

# c is linear contrast (could be filter response)

# from equation 6.5, p 205 of Felix's thesis:
crf <- function(c,param_list){
  eta<- param_list$eta
  kappa <- param_list$kappa
  beta <- param_list$beta
  alpha <- param_list$alpha
  r <- (alpha * c^eta) / (beta + c^kappa)
}


source(paste(getwd(),'/funs/calc_tvc.R',sep=""))


# set up params:
eta<- c(3, 4)
kappa <- c(2, 2.5)
beta <- c(5e-5,5e-6,5e-7)
alpha <- c(30)
alpha_ratio <- c(1) # dprime increment to test for tvc.

# simulate some curves fit from this function to examine the effects of various params.
contrasts <- c(0, exp(seq(log(0.001),log(0.2),length=200))) # e.g. in percent.
gen <- expand.grid(c = contrasts, eta= eta, kappa = kappa, beta = beta, alpha = alpha, alpha_ratio = alpha_ratio)
param_list <- list(eta= gen$eta, kappa = gen$kappa, beta = gen$beta, alpha = gen$alpha)
gen$r <- crf(c = gen$c, param_list = param_list)

fig <- ggplot(gen,aes(x=c,y=r,colour=factor(beta))) + facet_grid(eta~ kappa) + geom_line()
fig <- fig + scale_x_log10(name="Contrast")
fig <- fig + ylab("Response (d' units)")
fig <- fig + ggtitle('Contrast response function (CRF)')
fig <- fig + theme(legend.position="top")
fig
ggsave(paste0(getwd(),'/demos/wichmann_crf.pdf'),width=6,height=5)

#--------------------
# TvC functions:

gen$thresh <- calc_tvc(crf_fun=crf,c=gen$c, param_list=param_list, d_prime_inc=gen$alpha_ratio[1])

fig <- ggplot(gen,aes(x=c,y=thresh,colour=factor(beta))) + facet_grid(eta~ kappa) + geom_line()
fig <- fig + scale_x_log10(name="Contrast (%)")
fig <- fig + scale_y_log10("Threshold (% contrast)")
fig <- fig + ggtitle('Threshold versus Contrast (TvC)')
fig <- fig + ggtitle('Threshold versus contrast')
fig <- fig + theme(legend.position="top")
fig
ggsave(paste(getwd(),'/demos/wichmann_tvc.pdf',sep=""),width=6,height=5)




