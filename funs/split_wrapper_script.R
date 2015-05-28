## wrapper script to allow fast fitting using same settings.

rm(list=ls())


#-------------------
# set split factor:

# split_factor_string = 'em_prevSacc_direction_factor'
# split_factor_string = 'em_nextSacc_time_StimOnToOn_factor'
# split_factor_string = 'em_cumDist_factor'
# split_factor_string = 'H_2D_40_factor'
split_factor_string = 'H_3D_00_factor'

#-------------------
# set model name:

model_name = 'crf_4_param'
# model_name = 'csf_model'

#-------------------
# do fits, plots:

source(paste(getwd(),'/crf_fitting/split_wrapper_funs.R',sep=""))

fit_split(model_name = model_name,split_factor_string=split_factor_string,iter = 10000, thin = 20)
plot_split(model_name = model_name,split_factor_string=split_factor_string)
