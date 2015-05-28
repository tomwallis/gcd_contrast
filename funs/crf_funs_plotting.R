## functions to plot contrast response functions for mcmc estimates derived from Stan.
# TSAW, contact tsawallis@gmail.com.

# plotting functions will produce y, ymin and ymax values for a given fit object and x-values...


#------------------------
# generate a data frame containing predictions for different curves:

generate_predictions <- function(fit,data_frame,model_name,log_contrast_range='default',log_sf_range='default',mid_val='mean'){
  library(ggplot2)
  library(rstan)
  
  source('~/Dropbox/R_Functions/mcmc_helper_funcs.R')
  source(paste(getwd(),'/functions/calc_tvc.R',sep=""))
  source(paste(getwd(),'/functions/helpers.R',sep="")) 
  
  crf_fun <- crf_function(model_name)
  crf_sample <- sample_function(model_name)
  
  params <- extract(fit)
  n_samples <- length(params$lp__)
  
  # determine contrast over which to plot:
  if(any(log_contrast_range=='default')){
    min_contrast <- log10(min(data_frame$band_patch_Target[data_frame$band_patch_Target>0]))
    max_contrast <- log10(max(data_frame$band_patch_Target))
    contrasts <- c(0,10^seq(min_contrast,max_contrast,length=100))
  } else {
    contrasts <- c(0,10^seq(log_contrast_range[1],log_contrast_range[2],length=100))
  }
  
  # determine sf over which to plot csf:
  if(any(log_sf_range=='default')){
    min_sf <- log10(0.375)
    max_sf <- log10(24)
    sf_plot_x <- 10^seq(min_sf,max_sf,length=100)
  } else {
    sf_plot_x <- 10^seq(log_sf_range[1],log_sf_range[2],length=100)
  }
  
  alphas <- sort(unique(data_frame$alpha))
  
  # set up empty data frames in which to place predictions:
  pred_frame <- expand.grid(subject=levels(data_frame$subject),c=contrasts,alphas=alphas,sf_factor=levels(data_frame$sf_factor),d_prime_inc=c(1))
  
  pred_frame$crf <- NA
  pred_frame$crf_ymin <- NA
  pred_frame$crf_ymax <- NA
  pred_frame$tvc <- NA
  pred_frame$tvc_ymin <- NA
  pred_frame$tvc_ymax <- NA
  
  csf_frame <- expand.grid(subject=levels(data_frame$subject),sf=sf_plot_x)
  csf_frame$y <- NA
  csf_frame$ymin <- NA
  csf_frame$ymax <- NA    
  
  
  for (subj in 1:length(levels(data_frame$subject))){
    this_subj <- levels(data_frame$subject)[subj]
    
    if(model_name=='csf_model' | model_name=='csf_model_2'){
      #---------------- 
      # generate CSFs:
      # sf doesn't matter here since z not being used.
      param_list <- extract_params(params,model_name,subj_number=subj,sf_value=1) 
      
      ### csf curve for every mcmc sample:
      csf <- sapply(1:n_samples,yqm_csf_sample,f=sf_plot_x,param_list)
      csf <- t(csf)
      if(mid_val == 'mean') y <- colMeans(csf)
      if(mid_val == 'median') y <- vectorQuants(csf,probs=0.5)
      hdis <- vectorHDI(csf)
      
      csf_frame$y[csf_frame$subject==this_subj] <- y
      csf_frame$ymin[csf_frame$subject==this_subj] <- hdis[1,]
      csf_frame$ymax[csf_frame$subject==this_subj] <- hdis[2,]
    }
    
    for (sf in 1:length(levels(data_frame$sf_factor))){
      this_sf <- levels(data_frame$sf_factor)[sf]
      this_sf_num <- sort(unique(data_frame$sf))[sf]
      
      if(model_name=='csf_model' | model_name=='csf_model_2'){
        param_list <- extract_params(params,model_name,subj_number=subj,sf_value=this_sf_num)
      } else {
        param_list <- extract_params(params,model_name,subj_number=subj,sf_number=sf)
      }
      
      ### crf curve for every mcmc sample:
      crfs <- sapply(1:n_samples,crf_sample,c=contrasts,param_list)
      crfs <- t(crfs)
      if(mid_val == 'mean') y <- colMeans(crfs)
      if(mid_val == 'median') y <- vectorQuants(crfs,probs=0.5)
      hdis <- vectorHDI(crfs)
      
      pred_frame$crf[pred_frame$subject==this_subj & pred_frame$sf_factor==this_sf] <- y
      pred_frame$crf_ymin[pred_frame$subject==this_subj & pred_frame$sf_factor==this_sf] <- hdis[1,]
      pred_frame$crf_ymax[pred_frame$subject==this_subj & pred_frame$sf_factor==this_sf] <- hdis[2,]
      
      ### TvC function
      for (inc in 1:length(unique(pred_frame$d_prime_inc))){
        this_inc <- unique(pred_frame$d_prime_inc)[inc]
        thresh_samples <- sapply(1:n_samples,calc_tvc_vector,
                                 crf_sample=crf_sample,c=contrasts,param_list=param_list,
                                 d_prime_inc=this_inc)
        thresh_samples <- t(thresh_samples)
        if(mid_val == 'mean') y <- colMeans(thresh_samples)
        if(mid_val == 'median') y <- vectorQuants(thresh_samples,probs=0.5)
        hdis <- vectorHDI(thresh_samples)
        
        this_rows <- pred_frame$subject==this_subj & pred_frame$sf_factor==this_sf & pred_frame$d_prime_inc==this_inc
        
        pred_frame$tvc[this_rows] <- y
        pred_frame$tvc_ymin[this_rows] <- hdis[1,]
        pred_frame$tvc_ymax[this_rows] <- hdis[2,]
      }
    }
  }
  
  #-------------------
  # generate other structures to plot params.
  out <- generate_param_frames(params,data_frame,model_name)
  
  out_list <- list(pred_frame=pred_frame,csf_frame=csf_frame,thresh_frame=out$thresh_frame,param_frame=out$param_frame,csf_param_frame=out$csf_param_frame)
  return(out_list)
}


#---------------------
# plot CRFs at each spatial frequency:
plot_crfs <- function(generated_list,data_frame,output_file='none',log_contrast_range='default',plot_dims=c(6,6)){
  source('~/Dropbox/R_Functions/pleasingColours.R')
  library(scales)
  
  plot_frame <- generated_list$pred_frame
  
  # sample data_frame to reduce rug size.
  samp <- data_frame[sample(1:nrow(data_frame),size=2000),]
  
  # are we plotting with a split factor?
  if('split_factor' %in% colnames(plot_frame) == TRUE){
    fig <- ggplot(plot_frame,aes(x=c,y=crf,colour=split_factor,fill=split_factor))
    fig <- fig + facet_grid(subject ~ sf_factor)
    fig <- fig + geom_ribbon(aes(ymin=crf_ymin,ymax=crf_ymax,colour=NULL),alpha=.2)
    fig <- fig + geom_line()
    fig <- fig + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
    fig <- fig + xlab('Contrast [band energy]') + ylab("Response [d' units]")
    fig <- fig + theme_grey(base_size=11)
    
    # zoom to log_contrast_range:
    if(any(log_contrast_range!='default'))
      fig <- fig + coord_cartesian(xlim=10^log_contrast_range)
    
    # fix legend:
    fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
    fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
    fig <- fig + theme(legend.position='top')
  }
  
  if('split_factor' %in% colnames(plot_frame) == FALSE){
    fig <- ggplot(plot_frame,aes(x=c,y=crf))
    fig <- fig + facet_grid(subject ~ sf_factor)
    fig <- fig + geom_ribbon(aes(ymin=crf_ymin,ymax=crf_ymax,colour=NULL),alpha=.2)
    fig <- fig + geom_line()
    fig <- fig + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))   
    fig <- fig + geom_rug(data=samp,aes(x=band_patch_Target,y=NULL,colour=NULL),sides='b',alpha=.1)
    fig <- fig + xlab('Contrast [band energy]') + ylab("Response [d' units]")
    fig <- fig + theme_grey(base_size=11)
    
    # zoom to log_contrast_range:
    if(any(log_contrast_range!='default'))
      fig <- fig + coord_cartesian(xlim=10^log_contrast_range)
  }
  
  if(output_file!='none'){
    ggsave(file=paste(getwd(),'/crf_fitting/',output_file,sep=""),width=plot_dims[1],height=plot_dims[2])  
  } else {
    print(fig)
  }  
}

#---------------------
# plot all crf params:
plot_params <- function(generated_list,data_frame,output_file='none',plot_dims=c(3.5,5)){
  source('~/Dropbox/R_Functions/pleasingColours.R')
  library(scales)
  
  plot_frame <- generated_list$param_frame
  
  # are we plotting with a split factor?
  if('split_factor' %in% colnames(plot_frame) == TRUE){
    # plot other parameter credible distributions:
    fig <- ggplot(plot_frame,aes(x=value,fill=split_factor)) + geom_histogram(alpha=0.7,position='identity')
    fig <- fig + facet_grid(subject ~ param,scales="free_x")
    fig <- fig + xlab('Parameter Value') + ylab('Frequency')
    fig <- fig + theme_grey(base_size=11)
    
    # fix legend:
    fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
    fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
    fig <- fig + theme(legend.position='top')
  }
  
  if('split_factor' %in% colnames(plot_frame) == FALSE){
    # plot other parameter credible distributions:
    fig <- ggplot(plot_frame,aes(x=value)) + geom_histogram()
    fig <- fig + facet_grid(subject ~ param,scales="free_x")
    fig <- fig + xlab('Parameter Value') + ylab('Frequency')
    fig <- fig + theme_grey(base_size=11)
  }
  
  if(output_file!='none'){
    ggsave(file=paste(getwd(),'/crf_fitting/',output_file,sep=""),width=plot_dims[1],height=plot_dims[2])  
  } else {
    print(fig)
  }  
}

#---------------------
# plot thresholds:
plot_threshs <- function(generated_list,data_frame,output_file='none',plot_dims=c(4,4)){
  source('~/Dropbox/R_Functions/pleasingColours.R')
  source('~/Dropbox/R_Functions/mcmc_helper_funcs.R')
  library(plyr)
  library(scales)
  
  plot_frame <- generated_list$thresh_frame
  
  if('split_factor' %in% colnames(plot_frame) == TRUE){
    # do some pre-summaries to plot mean and hdi:
    sum_frame <- ddply(plot_frame,.(subject,sf_factor,split_factor),summarise,y=mean(thresh),ymin=hdiOfMCMC(thresh)[1],ymax=hdiOfMCMC(thresh)[2])
    
    fig <- ggplot(sum_frame,aes(x=sf_factor,y=1/y,colour=split_factor))
    fig <- fig + geom_linerange(aes(ymin=1/ymin,ymax=1/ymax),position=position_dodge(width=.3))
    fig <- fig + geom_point(position=position_dodge(width=.3))
    fig <- fig + facet_wrap(~ subject, ncol=3)
    fig <- fig + scale_y_log10(labels = trans_format("log10", math_format(10^.x)))
    fig <- fig + xlab('Spatial frequency band') + ylab('Sensitivity [1 / z]')
    fig <- fig + theme_grey(base_size=11)
    fig <- fig + theme(axis.text.x = element_text(angle=60,vjust=.7))
    
    # fix legend:
    fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
    fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
    fig <- fig + theme(legend.position='top')
  }
  
  if('split_factor' %in% colnames(plot_frame) == FALSE){
    # do some pre-summaries to plot mean and hdi:
    sum_frame <- ddply(plot_frame,.(subject,sf_factor),summarise,y=mean(thresh),ymin=hdiOfMCMC(thresh)[1],ymax=hdiOfMCMC(thresh)[2])
    
    fig <- ggplot(sum_frame,aes(x=sf_factor,y=1/y))
    fig <- fig + geom_linerange(aes(ymin=1/ymin,ymax=1/ymax))
    fig <- fig + geom_point()
    fig <- fig + facet_wrap(~ subject, ncol=3)
    fig <- fig + scale_y_log10(labels = trans_format("log10", math_format(10^.x)))
    fig <- fig + xlab('Spatial frequency band') + ylab('Contrast sensitivity (1/threshold)')
    fig <- fig + theme_grey(base_size=11)
    fig <- fig + theme(axis.text.x = element_text(angle=60,vjust=.7))
  }
  
  if(output_file!='none'){
    ggsave(file=paste(getwd(),'/crf_fitting/',output_file,sep=""),width=plot_dims[1],height=plot_dims[2])  
  } else {
    print(fig)
  }  
}

#---------------------
# plot TvCs at each spatial frequency:
plot_tvcs <- function(generated_list,data_frame,output_file='none',log_contrast_range='default',plot_dims=c(6,6)){
  source('~/Dropbox/R_Functions/pleasingColours.R')
  library(scales)
  
  plot_frame <- generated_list$pred_frame
  
  plot_frame$d_prime_inc <- factor(plot_frame$d_prime_inc)
  
  # sample data_frame to reduce rug size.
  samp <- data_frame[sample(1:nrow(data_frame),size=2000),]
  
  if('split_factor' %in% colnames(plot_frame) == TRUE){
    if (length(levels(plot_frame$d_prime_inc))==1)
      fig <- ggplot(plot_frame,aes(x=c,y=tvc,colour=split_factor,fill=split_factor))  
    if (length(levels(plot_frame$d_prime_inc))>1)
      fig <- ggplot(plot_frame,aes(x=c,y=tvc,linetype=d_prime_inc,colour=split_factor,fill=split_factor))  
    fig <- fig + facet_grid(subject ~ sf_factor)
    fig <- fig + geom_ribbon(aes(ymin=tvc_ymin,ymax=tvc_ymax,colour=NULL),alpha=.2)
    fig <- fig + geom_line()
    fig <- fig + geom_rug(data=samp,aes(x=band_patch_Target,y=NULL,colour=NULL,fill=NULL),sides='b',alpha=.1)
    fig <- fig + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
    fig <- fig + scale_y_log10()
    fig <- fig + theme_grey(base_size=11)
    fig <- fig + theme(legend.position='top')
    fig <- fig + scale_linetype_discrete(name="d' increment")
    fig <- fig + xlab('Contrast [band energy]') + ylab("Threshold")
    
    # fix legend:
    fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
    fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
    fig <- fig + theme(legend.position='top')
    
    # zoom to log_contrast_range:
    if(any(log_contrast_range!='default'))
      fig <- fig + coord_cartesian(xlim=10^log_contrast_range)
  }
  
  if('split_factor' %in% colnames(plot_frame) == FALSE){
    if (length(levels(plot_frame$d_prime_inc))==1)
      fig <- ggplot(plot_frame,aes(x=c,y=tvc))  
    if (length(levels(plot_frame$d_prime_inc))>1)
      fig <- ggplot(plot_frame,aes(x=c,y=tvc,linetype=d_prime_inc))  
    fig <- fig + facet_grid(subject ~ sf_factor)
    fig <- fig + geom_ribbon(aes(ymin=tvc_ymin,ymax=tvc_ymax,colour=NULL),alpha=.2)
    fig <- fig + geom_line()
    fig <- fig + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
    fig <- fig + scale_y_log10()
    fig <- fig + geom_rug(data=samp,aes(x=band_patch_Target,y=NULL,colour=NULL,linetype=NULL),sides='b',alpha=.1)
    fig <- fig + theme_grey(base_size=11)
    fig <- fig + theme(legend.position='top')
    fig <- fig + scale_linetype_discrete(name="d' increment")
    fig <- fig + xlab('Contrast [band energy]') + ylab("Threshold")
    
    # zoom to log_contrast_range:
    if(any(log_contrast_range!='default'))
      fig <- fig + coord_cartesian(xlim=10^log_contrast_range)
  }
  
  if(output_file!='none'){
    ggsave(file=paste(getwd(),'/crf_fitting/',output_file,sep=""),width=plot_dims[1],height=plot_dims[2])  
  } else {
    print(fig)
  }  
}

#---------------------
# plot csf parameters:
plot_csf_params <- function(generated_list,data_frame,output_file='none',plot_dims=c(5.5,4)){
  source('~/Dropbox/R_Functions/pleasingColours.R')
  library(plyr)
  library(scales)
  
  plot_frame <- generated_list$csf_param_frame
  
  # are we plotting with a split factor?
  if('split_factor' %in% colnames(plot_frame) == TRUE){
    # plot other parameter credible distributions:
    fig <- ggplot(plot_frame,aes(x=value,fill=split_factor)) + geom_histogram(alpha=0.7,position='identity')
    fig <- fig + facet_grid(subject ~ param,scales="free_x")
    fig <- fig + xlab('Parameter Value') + ylab('Frequency')
    fig <- fig + theme_grey(base_size=11)
    
    # fix legend:
    fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
    fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
    fig <- fig + theme(legend.position='top')
    fig <- fig + theme_grey(base_size=11)
  }
  
  if('split_factor' %in% colnames(plot_frame) == FALSE){
    # plot other parameter credible distributions:
    fig <- ggplot(plot_frame,aes(x=value)) + geom_histogram()
    fig <- fig + facet_grid(subject ~ param,scales="free_x")
    fig <- fig + xlab('Parameter Value') + ylab('Frequency')
    fig <- fig + theme_grey(base_size=11)
  }
  
  if(output_file!='none'){
    ggsave(file=paste(getwd(),'/crf_fitting/',output_file,sep=""),width=plot_dims[1],height=plot_dims[2])  
  } else {
    print(fig)
  }  
}

#---------------------
# plot sensitivity at each spatial frequency:
plot_csf <- function(generated_list,data_frame,output_file='none',log_sf_range='default',plot_dims=c(2.5,5.5)){
  source('~/Dropbox/R_Functions/pleasingColours.R')
  library(scales)
  
  plot_frame <- generated_list$csf_frame
  
  breaks <- c(0.5, 1, 2, 4, 8, 16)
  
  # are we plotting with a split factor?
  if('split_factor' %in% colnames(plot_frame) == TRUE){
    fig <- ggplot(plot_frame,aes(x=sf,y=y,colour=split_factor,fill=split_factor))
    fig <- fig + facet_wrap(~ subject,ncol=1)
    fig <- fig + geom_ribbon(aes(ymin=ymin,ymax=ymax,colour=NULL),alpha=.2)
    fig <- fig + geom_line()
    fig <- fig + scale_x_log10(breaks=breaks)
    fig <- fig + xlab('Spatial Frequency [c / deg]') + ylab("Sensitivity")
    fig <- fig + theme_grey(base_size=11)
    
    # zoom to log_contrast_range:
    if(any(log_sf_range!='default'))
      fig <- fig + coord_cartesian(xlim=10^log_sf_range)
    
    # fix legend:
    fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
    fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
    fig <- fig + theme(legend.position='top')
  }
  
  if('split_factor' %in% colnames(plot_frame) == FALSE){
    fig <- ggplot(plot_frame,aes(x=sf,y=y))
    fig <- fig + facet_wrap(~ subject,ncol=1)
    fig <- fig + geom_ribbon(aes(ymin=ymin,ymax=ymax,colour=NULL),alpha=.2)
    fig <- fig + geom_line()
    fig <- fig + scale_x_log10(breaks=breaks)
    fig <- fig + xlab('Spatial Frequency [c / deg]') + ylab("Sensitivity")
    fig <- fig + theme_grey(base_size=11)
    
    # zoom to log_contrast_range:
    if(any(log_sf_range!='default'))
      fig <- fig + coord_cartesian(xlim=10^log_sf_range)
  }
  
  if(output_file!='none'){
    ggsave(file=paste(getwd(),'/crf_fitting/',output_file,sep=""),width=plot_dims[1],height=plot_dims[2])  
  } else {
    print(fig)
  }  
}

