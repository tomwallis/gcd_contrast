## produce plots to check MCMC convergence etc.

rm(list=ls())

library(ggmcmc)

# GLM fits ---------------------------------------

load(paste0(getwd(),'/output/full_glm.RData'))
print(fit)

ggmcmc(ggs(fit), file = 'figs/full_glm_traceplots.pdf', plot = c('ggs_density()','ggs_autocorrelation()'))
