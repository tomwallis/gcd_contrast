## wrapper script to allow fast fitting using same settings.

rm(list=ls())

# model_name = 'crf_2_param'
# model_name = 'crf_3_param'
# model_name = 'crf_4_param'
# model_name = 'crf_full_model'
# model_name = 'csf_model'
model_name = 'csf_model_2'


source(paste(getwd(),'/crf_fitting/wrapper_funs.R',sep=""))

fit_wrapper(model_name = model_name,iter = 10000, thin = 20)
plot_wrapper(model_name = model_name)
