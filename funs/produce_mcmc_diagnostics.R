## produce plots to check MCMC convergence etc.

rm(list=ls())

library(ggmcmc)

# GLM fits ---------------------------------------

load(paste0(getwd(),'/output/single_level_glm.RData'))
simple_fit <- fit
print(fit)

load(paste0(getwd(),'/output/multilevel_glm.RData'))
multilevel_fit <- fit
print(fit)

rm(fit)

ggmcmc(ggs(simple_fit), file = 'figs/single_level_glm_traceplots.pdf', plot = c('ggs_density()','ggs_autocorrelation()'))

ggmcmc(ggs(multilevel_fit), file = 'figs/multilevel_glm_traceplots.pdf', plot = c('ggs_density()','ggs_autocorrelation()'))


# transducer fits ---------------------------------------

load(paste0(getwd(),'/output/single_level_transducer.RData'))
transducer_a <- fit
print(fit)

load(paste0(getwd(),'/output/single_level_transducer_B.RData'))
transducer_b <- fit
print(fit)

rm(fit)

ggmcmc(ggs(transducer_a), file = 'figs/single_level_transducer_traceplots.pdf', plot = c('ggs_density()','ggs_autocorrelation()'))

ggmcmc(ggs(transducer_b), file = 'figs/single_level_transducer_B_traceplots.pdf', plot = c('ggs_density()','ggs_autocorrelation()'))
