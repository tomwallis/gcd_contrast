## plot crossvalidation results against one another.
# TSAW.

rm(list=ls())

library(ggplot2)
library(plyr)
library(reshape2)
source(paste(getwd(),'/functions/cross_validation.R',sep=""))
source('~/Dropbox/R_Functions/pleasingColours.R')

model_filename <- c('crf_2_param','crf_4_param','csf_model','GLM_')
model_name <- c('p & q fixed','all free','CSF Model','GLM')

curve_frame <- data.frame()
label_frame <- data.frame()
y_hat_frame <- data.frame()
aucs <- data.frame()
devs <- data.frame()

#-----------------
# load roc output for each model, place into a data frame:
for (i in 1 : length(model_filename)){
  load(paste(getwd(),'/crf_crossval/',model_filename[i],'cv_fit.RData',sep=""))
  this_frame <- data.frame(model=model_name[i],x=roc$x,y=roc$y,auc=roc$auc)
  curve_frame <- rbind(curve_frame,this_frame)
  
#   # cumulative matrix of model predictions:
  this_frame <- data.frame(model=model_name[i],y_hat=y_hat,sample_n=1:(length(y_hat)))
  y_hat_frame <- rbind(y_hat_frame,this_frame)
  
  this_frame <- data.frame(model=model_name[i],x = 0.65, y = NA, label = paste(model_name[i],', AUC = ',roc$auc,sep=""), auc = roc$auc)
  label_frame <- rbind(label_frame,this_frame)
  
  # cumulative data frame of area under curve by subject:
  aucs <- rbind(aucs,auc_frame)  
  devs <- rbind(devs,dev_frame)
}

#-----------------
# sort labels by order of auc, assign y position of label:
curve_frame$model <- reorder(curve_frame$model, -curve_frame$auc)
label_frame$model <- reorder(label_frame$model, -label_frame$auc)
aucs$model <- reorder(aucs$model, -aucs$auc)
devs$model <- reorder(devs$model, -devs$dev)
label_frame <- label_frame[order(label_frame$auc),]

for (i in 1:length(levels(label_frame$model))){
  label_frame$y[i] <- 0.1 + (0.05*i)
}

# order factor with GLM at the end:
curve_frame$model <- factor(curve_frame$model,levels = c('all free','CSF model', 'p & q fixed','GLM'))
aucs$model <- factor(aucs$model,levels = c('all free','CSF Model', 'p & q fixed','GLM'))
devs$model <- factor(devs$model,levels = c('all free','CSF Model', 'p & q fixed','GLM'))

# combine auc and devs to plot together.
combined_frame <- aucs
combined_frame$measure <- 'AUC'
combined_frame <- rename(combined_frame, replace=c('auc' = 'y'))
devs$measure <- 'Deviance'
devs <- rename(devs, replace=c('dev' = 'y'))
combined_frame <- rbind(combined_frame,devs)
#-----------------
# plot results:
fig <- ggplot(curve_frame,aes(x=x,y=y,colour=model)) + geom_line()
fig <- fig + geom_line(aes(x=x,y=x),linetype='dashed',colour="black")
fig <- fig + xlab('False positive rate') + ylab('True positive rate')
fig <- fig + geom_text(data=label_frame,aes(x=x,y=y,label=label),size=3)
fig <- fig + scale_colour_manual(values=cbbPalette[-1],guide=FALSE)
fig <- fig + theme_grey(base_size=11)
# fig
ggsave(file=paste(getwd(),'/crf_crossval/cv_results.pdf',sep=""),width=3,height=3)

#-----------------
# plot auc by subject, model:

fig <- ggplot(aucs,aes(x=subject,y=auc,colour=model,fill=model))
fig <- fig + geom_bar(stat="identity",position=position_dodge())
fig <- fig + xlab('Subject') + ylab('Area under ROC')
fig <- fig + scale_colour_manual(name='',values=cbbPalette[-1])
fig <- fig + scale_fill_manual(name='',values=cbbPalette[-1])
fig <- fig + coord_cartesian(ylim=c(0.58,0.74))
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
# fig
ggsave(file=paste(getwd(),'/crf_crossval/cv_barplots.pdf',sep=""),width=4,height=3.5)

#-----------------
# plot deviance by subject, model

# only full data deviance:
devs <- subset(devs,subject=='All data')

fig <- ggplot(devs,aes(x=model,y=y))
fig <- fig + geom_bar(stat='identity')
fig <- fig + xlab('Model') + ylab('Crossvalidated deviance')
fig <- fig + theme_grey(base_size=11)
fig <- fig + theme(legend.position='top')
fig <- fig + coord_cartesian(ylim=c(min(devs$y)*0.98,max(devs$y)*1.02))
fig <- fig + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# fig
ggsave(file=paste(getwd(),'/crf_crossval/cv_deviance.pdf',sep=""),width=4,height=3.5)
