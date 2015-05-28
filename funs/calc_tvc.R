## Method to find a numerical appoximation for a TvC function. 
# Based on Andrew Haun's matlab code; adapted by TSAW.

calc_tvc <- function(crf_fun,c=10^seq(-3,0,length=100),param_list=list(p=2,q=.4,z=.01,rmax=30),d_prime_inc=1,grad_step=.01){
  
  r_standard <- crf_fun(c,param_list)
  r_increment <- crf_fun((c + grad_step),param_list)
  
  thresh <- d_prime_inc * grad_step / (r_increment - r_standard)
  return(thresh)
}

#---------------------
# same, with a sample index for use with sapply.
calc_tvc_vector <- function(i,crf_sample,c=10^seq(-3,0,length=100),
                            param_list=list(p=c(1.5,2),q=c(.2,.4),z=c(.01,.1),rmax=c(5,30)),
                            d_prime_inc=1,grad_step=.01){
  # vectorised, to be used with sapply.
  r_standard <- crf_sample(i,c,param_list)
  r_increment <- crf_sample(i,(c + grad_step),param_list)
  thresh <- d_prime_inc * grad_step / (r_increment - r_standard)
  return(thresh)
}
