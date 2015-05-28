## wrapper functions for fitting split data plots.

#------------------------
fit_split <- function(model_name,split_factor_string,iter=1000,thin=5){
  # load functions:
  source(paste(getwd(),'/functions/cross_validation.R',sep=""))
  source(paste(getwd(),'/functions/crf_funs_fitting.R',sep=""))
  source(paste(getwd(),'/functions/helpers.R',sep=""))
  
  # load data frame:
  load(paste(getwd(),'/crf_fitting/split_factor_results/csf_data_split_factors.RData',sep=''))
  
  # process to xy:
  xy <- process_df_to_xy(dat, model_name = model_name, na_omit = TRUE,split_factor=split_factor_string)
  x <- xy$x
  y <- xy$y
  rm(xy)
  
  # how many levels does the split factor have?
  text <- paste('factor_levels <- levels(x$',split_factor_string,')',sep="")
  eval(parse(text=text))
  
  # create a grouping variable, separate data set for each level, fit:
  for (g in 1 : length(factor_levels)){
    text <- paste('group',g,' <- x$',split_factor_string,'=="',factor_levels[g],'"',sep="")
    eval(parse(text=text))
    
    text <- paste('x',g,' <- x[group',g,'==TRUE,]',sep="")
    eval(parse(text=text))
    
    text <- paste('y',g,' <- y[group',g,'==TRUE]',sep="")
    eval(parse(text=text))
    
    # fit models:
    text <- paste('fit',g,' <- fit_crf_parallel(x',g,',y',g,', iter=iter, model_name = model_name, thin = thin, n_chains = 4)',sep="")
    eval(parse(text=text))
    
    if(g == 1){
      fit_string <- 'fit1'
    } else {
      fit_string <- paste(fit_string,',fit',g,sep="")
    }
  }
  
  # save fits:
  text <- paste('save(',fit_string,',file="',getwd(),'/crf_fitting/split_factor_results/',split_factor_string,'_',model_name,'_fit.RData")',sep="")
  eval(parse(text=text))
  
}


#------------------------
plot_split <- function(model_name,split_factor_string){
  # load functions:
  source(paste(getwd(),'/functions/cross_validation.R',sep=""))
  source(paste(getwd(),'/functions/crf_funs_plotting.R',sep=""))
  source(paste(getwd(),'/functions/helpers.R',sep=""))

  # load fits:
  load(paste(getwd(),'/crf_fitting/split_factor_results/',split_factor_string,'_',model_name,'_fit.RData',sep=''))
  
  # load data frame:
  load(paste(getwd(),'/crf_fitting/split_factor_results/csf_data_split_factors.RData',sep=''))
  
  # process nas out of data frame:
  xy <- process_df_to_xy(dat, model_name = model_name, na_omit = TRUE,split_factor=split_factor_string)
  x <- xy$x
  y <- xy$y
  rm(xy)
  dat <- cbind(x,y)
  
  # set up for plotting:
  log_cont_range <- c(-3.2,-0.5)
  log_sf_range <- log10(c(0.375,24))
  
  # what are the split factor levels?
  text <- paste('factor_levels <- levels(x$',split_factor_string,')',sep="")
  eval(parse(text=text))
  
  # create a grouping variable, separate data set for each level, generate predictions separately:
  for (g in 1 : length(factor_levels)){
    text <- paste('group',g,' <- x$',split_factor_string,'=="',factor_levels[g],'"',sep="")
    eval(parse(text=text))
    
    text <- paste('x',g,' <- x[group',g,'==TRUE,]',sep="")
    eval(parse(text=text))
    
    text <- paste('y',g,' <- y[group',g,'==TRUE]',sep="")
    eval(parse(text=text))
    
    # bind back into data set:
    text <- paste('dat',g,' <- cbind(x',g,',y',g,')',sep="")
    eval(parse(text=text))
    
    # generate predictions:
    text <- paste('out_list',g,' <- generate_predictions(fit',g,',dat',g,',model_name="',model_name,'",log_contrast_range=c(',log_cont_range[1],',',log_cont_range[2],'),log_sf_range=c(',log_sf_range[1],',',log_sf_range[2],'))',sep="")
    eval(parse(text=text))
}
  
  # separate predictions must be combined into one data frame for plotting. out_lists contain a number of sub-frames that are required by the subsequent plotting functions. Do each one.
  for (g in 1 : length(factor_levels)){
    # add split factor field to each sub-frame in this out_list:
    text <- paste('out_list',g,'$pred_frame$split_factor <- as.factor(factor_levels[',g,'])',sep="")
    eval(parse(text=text))
    
    text <- paste('out_list',g,'$csf_frame$split_factor <- as.factor(factor_levels[',g,'])',sep="")
    eval(parse(text=text))
    
    text <- paste('out_list',g,'$thresh_frame$split_factor <- as.factor(factor_levels[',g,'])',sep="")
    eval(parse(text=text))
    
    text <- paste('out_list',g,'$param_frame$split_factor <- as.factor(factor_levels[',g,'])',sep="")
    eval(parse(text=text))
    
    if(length(out_list1$csf_param_frame)!=0){
      text <- paste('out_list',g,'$csf_param_frame$split_factor <- as.factor(factor_levels[',g,'])',sep="")
      eval(parse(text=text))
      }
  }
  
  
  # bind all the frames together:
  out_list <- list()
  string1 <- 'out_list$pred_frame <- rbind(out_list1$pred_frame,'
  string2 <- 'out_list$csf_frame <- rbind(out_list1$csf_frame,'
  string3 <- 'out_list$thresh_frame <- rbind(out_list1$thresh_frame,'
  string4 <- 'out_list$param_frame <- rbind(out_list1$param_frame,'
  string5 <- 'out_list$csf_param_frame <- rbind(out_list1$csf_param_frame,'

  for (g in 2 : length(factor_levels)){
    if (g==length(factor_levels)){
      string1 <- paste(string1,'out_list',g,'$pred_frame)',sep="")
      string2 <- paste(string2,'out_list',g,'$csf_frame)',sep="")
      string3 <- paste(string3,'out_list',g,'$thresh_frame)',sep="")
      string4 <- paste(string4,'out_list',g,'$param_frame)',sep="")
      string5 <- paste(string5,'out_list',g,'$csf_param_frame)',sep="")
    }    
    else {
      string1 <- paste(string1,'out_list',g,'$pred_frame,',sep="")
      string2 <- paste(string2,'out_list',g,'$csf_frame,',sep="")
      string3 <- paste(string3,'out_list',g,'$thresh_frame,',sep="")
      string4 <- paste(string4,'out_list',g,'$param_frame,',sep="")
      string5 <- paste(string5,'out_list',g,'$csf_param_frame,',sep="")
    }
  }

  eval(parse(text=string1))
  eval(parse(text=string2))
  eval(parse(text=string3))
  eval(parse(text=string4))
  if(length(out_list1$csf_param_frame)!=0) eval(parse(text=string5))
  
  #-----------------
  # run plot commands. split factors should be detected automatically.
  plot_crfs(out_list,dat,output_file = paste('/split_factor_results/',split_factor_string,'_',model_name,'_crf.pdf',sep=""),log_contrast_range=log_cont_range)
  plot_tvcs(out_list,dat,output_file = paste('/split_factor_results/',split_factor_string,'_',model_name,'_tvc.pdf',sep=""),log_contrast_range=log_cont_range)
  plot_params(out_list,dat,output_file= paste('/split_factor_results/',split_factor_string,'_',model_name,'_params.pdf',sep=""))
  plot_threshs(out_list,dat,output_file= paste('/split_factor_results/',split_factor_string,'_',model_name,'_threshs.pdf',sep=""))  
  
  if(model_name=='csf_model' | model_name=='csf_model_2'){
    plot_csf(out_list,dat,output_file= paste('/split_factor_results/',split_factor_string,'_',model_name,'_csf.pdf',sep=""))  
    plot_csf_params(out_list,dat,output_file= paste('/split_factor_results/',split_factor_string,'_',model_name,'_csf_params.pdf',sep=""))  
  }
}
