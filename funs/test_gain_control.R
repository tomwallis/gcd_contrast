## testing activity of gain control terms....

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
  rmax <- param_list$rmax
  z <- param_list$z
  
  denom <- c^p + z^p
  
  r <- rmax * (c^(p+q) / denom )
} 


crf_gc <- function(c, norm_pool, param_list){
  p <- param_list$p
  q <- param_list$q
  rmax <- param_list$rmax
  z <- param_list$z
  w_1 <- param_list$w_1 
  w_2 <- param_list$w_2 
  w_3 <- param_list$w_3 
  
  weighted <- NA
  
  for (i in 1 : 3){
    text <- paste0('w <- w_',i)
    eval(parse(text=text))
    weighted[,i] <- w * norm_pool[,i]^p
  }
  
  denom <- rowsum(weighted) + z^p
  
  r <- rmax * (c^(p+q) / denom )
}

# Test some typical values ------------------
# set up params:
p <- c(2)
q <- c(.4)
z <- c(.02)
rmax <- c(50)
alpha_ratio <- c(1) # dprime increment to test for tvc.

w_1 <- c(2, 1)
w_2 <- c(1)
w_3 <- c(0.5, 0.25)

# simulate some curves fit from this function to examine the effects of various params.
contrasts <- 10^seq(-3,0,length=200) # e.g. in percent.
gen <- expand.grid(c = contrasts, p = p, q = q, z = z, rmax = rmax, alpha_ratio = alpha_ratio,
                   w_1 = w_1, w_2 = w_2, w_3 = w_3,
                   model = c('normal','gc'))

# gain pool contrast (e.g. off filter response) is noisy contrast:
param_list <- list(p = gen$p, q = gen$q, z = gen$z, rmax = gen$rmax, 
                   w_1 = gen$w_1, w_2 = gen$w_2, w_3 = gen$w_3)
gen$r <- crf(c = gen$c, param_list = param_list)

norm_pool <- matrix(cbind(contrasts + rnorm(length(contrasts),sd=sd(contrasts)),
                          contrasts + rnorm(length(contrasts),sd=sd(contrasts)),
                          contrasts + rnorm(length(contrasts),sd=sd(contrasts))), ncol=3)

gc_r <- crf_gc(c = gen$c[gen$model=='gc'], 
               norm_pool = norm_pool,
               param_list = param_list)

gen$r[gen$model=='gc'] <- crf_gc(c = gen$c, 
                                 norm_pool = norm_pool,
                                 param_list = param_list)

break_locs <- c(.001,.002,.004,.01,.02,.04,.1,.2,.4,1)

fig <- ggplot(gen,aes(x=c,y=r,colour=model)) + facet_grid(w_1 ~ w_3) + geom_line()
fig <- fig + scale_x_log10(name="Contrast (%)",breaks=break_locs,labels=break_locs*100)
fig <- fig + ylab("Response (d' units)")
fig <- fig + ggtitle('Contrast response function (CRF)')
fig <- fig + theme(legend.position="top")
fig

