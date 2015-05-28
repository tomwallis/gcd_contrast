## wrapper functions to fit and plot a given model.
# TSAW

#------------------------
fit_wrapper <- function(model_name,iter=1000,thin=5){
  # load functions:
  source(paste(getwd(),'/functions/cross_validation.R',sep=""))
  source(paste(getwd(),'/functions/crf_funs_fitting.R',sep=""))
  source(paste(getwd(),'/functions/helpers.R',sep=""))
  
  # load data frame:
  load(paste(getwd(),'/data/csf_data_NAs.RData',sep=''))
  
  # process to xy:
  xy <- process_df_to_xy(dat, model_name = model_name, na_omit = TRUE)
  x <- xy$x
  y <- xy$y
  rm(xy)
  
  # fit model:
  fit <- fit_crf_parallel(x,y, iter=iter, model_name = model_name, thin = thin, n_chains = 4)
  
  # save fit and plot separately
  save(fit,file=paste(getwd(),'/crf_fitting/',model_name,'_fit.RData',sep=""))
}

#------------------------
plot_wrapper <- function(model_name){
  # load functions:
  source(paste(getwd(),'/functions/cross_validation.R',sep=""))
  source(paste(getwd(),'/functions/crf_funs_plotting.R',sep=""))
  source(paste(getwd(),'/functions/helpers.R',sep=""))
  
  # load fit object:
  load(paste(getwd(),'/crf_fitting/',model_name,'_fit.RData',sep=""))
  
  # load data frame:
  load(paste(getwd(),'/data/csf_data_NAs.RData',sep=''))
  
  # process nas out of data frame:
  xy <- process_df_to_xy(dat, model_name = model_name, na_omit = TRUE)
  x <- xy$x
  y <- xy$y
  rm(xy)
  dat <- cbind(x,y)
  
  log_cont_range <- c(-3.2,-0.5)
  log_sf_range <- log10(c(0.375,24))
  
  out_list <- generate_predictions(fit,dat,model_name=model_name,log_contrast_range=log_cont_range,log_sf_range=log_sf_range)
  plot_crfs(out_list,dat,output_file = paste(model_name,'_crf.pdf',sep=""),log_contrast_range=log_cont_range)
  plot_tvcs(out_list,dat,output_file = paste(model_name,'_tvc.pdf',sep=""),log_contrast_range=log_cont_range)
  plot_params(out_list,dat,output_file= paste(model_name,'_params.pdf',sep=""))
  plot_threshs(out_list,dat,output_file= paste(model_name,'_threshs.pdf',sep=""))  
  
  if(model_name=='csf_model' | model_name=='csf_model_2'){
    plot_csf(out_list,dat,output_file= paste(model_name,'_csf.pdf',sep=""))  
    plot_csf_params(out_list,dat,output_file= paste(model_name,'_csf_params.pdf',sep=""))  
  }
}

