## script to plot crf
# TSAW

rm(list=ls())

# load functions:
source(paste(getwd(),'/functions/crf_funs_plotting.R',sep=""))
source(paste(getwd(),'/functions/helpers.R',sep=""))

# load data frame:
load(paste(getwd(),'/data/csf_data_NAs.RData',sep=''))
# process nas out of data frame:
xy <- process_df_to_xy(dat, na_omit = TRUE)
x <- xy$x
y <- xy$y
rm(xy)
dat <- cbind(x,y)

#-----------------
# load fit objects in a loop, generate predictions:
model_names <- c('crf_2_param','crf_4_param','csf_model')
model_labels <- c('p & q fixed','all free','CSF model')

pred_frame <- data.frame()
thresh_frame <- data.frame()
param_frame <- data.frame()

for (i in 1 : length(model_names)){
  load(paste(getwd(),'/crf_fitting/',model_names[i],'_fit.RData',sep=''))  
  log_cont_range <- c(-3.2,-0.5) #'default' 
  out_list <- generate_predictions(fit,dat,model_name=model_names[i],log_contrast_range=log_cont_range)
  
  # squish predictions together for plotting:
  this_part <- out_list$pred_frame
  this_part$split_factor <- model_labels[i]
  pred_frame <- rbind(pred_frame,this_part)
  
  this_part <- out_list$thresh_frame
  this_part$split_factor <- model_labels[i]
  thresh_frame <- rbind(thresh_frame,this_part)
  
  this_part <- out_list$param_frame
  this_part$split_factor <- model_labels[i]
  param_frame <- rbind(param_frame,this_part)
}

pred_frame$split_factor <- factor(pred_frame$split_factor)
thresh_frame$split_factor <- factor(thresh_frame$split_factor)
param_frame$split_factor <- factor(param_frame$split_factor)

# squeeze back into out_list:
out_list <- list(pred_frame=pred_frame,thresh_frame=thresh_frame,param_frame=param_frame)

#-----------------
# plot using plot functions:
plot_crfs(out_list,dat,output_file = 'multi_crf.pdf',log_contrast_range=log_cont_range)
plot_tvcs(out_list,dat,output_file = 'multi_tvc.pdf',log_contrast_range=log_cont_range)
plot_params(out_list,dat,output_file= 'multi_params.pdf')
plot_threshs(out_list,dat,output_file= 'multi_threshs.pdf')
