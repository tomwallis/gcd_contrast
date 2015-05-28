# plot image statstic-related things -----------

rm(list = ls())

library(ggplot2)
library(mgcv)
library(plyr)
library(wutils)
library(psybayes)
library(reshape2)

# needs "full" data file.
load(paste0(getwd(),'/output/csf_data_full.RData'))
source(paste0(getwd(),"/funs/helpers.R"))

# exclude TW, AJ and LB (few trials)
dat <- dat[dat$subject!='AJ',]
dat <- dat[dat$subject!='LB',]
dat <- dat[dat$subject!='TSAW',]
dat$subject <- factor(dat$subject) # get rid of empty levels.

# rename subjects to anonymous codes:
levels(dat$subject) <- c('S1','S2','S3','S4','S5') # corresponds to MACD, PJB, LAL, CPT and AM.



# "Edge density" for several features at different scales, 3D features (i.e. spatiotemporal) --------------------- 


h <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
                            H_3D_00_target:H_3D_52_target))
h$invariant <- factor('H')
h <- melt(h, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'invariant'), variable.name = 'feature', value.name = 'value' )

s <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
                            S_3D_00_target:S_3D_52_target))
s$invariant <- factor('S')
s <- melt(s, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'invariant'), variable.name = 'feature', value.name = 'value' )
                            
k <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
                            K_3D_00_target:K_3D_52_target))
k$invariant <- factor('K')
k <- melt(k, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'invariant'), variable.name = 'feature', value.name = 'value' )

hsk <- rbind(h, s, k)
hsk <- na.omit(hsk)

# remove invariant name at the front of the "feature" string so that we can compare them all at the same scale:
remove_strings <- function(i){
  string <- substring(i, 3)
  return(string)
}

# unfortunately this is super slow:
hsk$feature <- sapply(hsk$feature, remove_strings)
hsk$feature <- factor(hsk$feature)

# take log of value, correcting for zeros:
hsk$log_value <- log_convert(hsk$value)

fig <- ggplot(hsk, aes(x = log_value, y = correct, colour=invariant, fill=invariant)) + 
  facet_wrap(~ feature, ncol = 3) + 
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5)) + 
  scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
#   geom_rug(sides = "b", alpha = 0.1) + 
  xlab('Feature value') + 
  scale_color_brewer(name='', palette='Dark2') +
  scale_fill_brewer(name='', palette='Dark2') +
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  theme(legend.position='top')
# fig
ggsave(file=paste0(getwd(),'/figs/image_features_full_3D.pdf'),width=5,height=8)




# based on checking out the full grid, just use K and 2 spatial scales at middle temporal scale:--------------------
sub_dat <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
                            K_3D_00_target:K_3D_52_target))
sub_dat <- melt(sub_dat, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct'), variable.name = 'feature', value.name = 'value' )
sub_dat <- na.omit(sub_dat)
sub_dat <- subset(sub_dat, feature == 'K_3D_51_target' | feature == 'K_3D_31_target' | feature == 'K_3D_01_target')
sub_dat$feature <- factor(sub_dat$feature)
sub_dat$feature <- factor(sub_dat$feature, levels = levels(sub_dat$feature)[c(3,2,1)])
levels(sub_dat$feature) <- c(".375-.75", "1.5-3", "12-24")

# take log of value, correcting for zeros:
sub_dat$log_val <- log_convert(sub_dat$value)

fig <- ggplot(sub_dat, aes(x = log_val, y = correct)) + 
  facet_wrap(~ feature, ncol = 3) + 
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=3),aes(fill=NULL), colour = 'black') + 
  geom_rug(sides = "b", alpha = 0.1) + 
  scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
  xlab('Feature value') + 
  theme_minimal(base_size=11) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
  theme(legend.position='top')
# fig
ggsave(file=paste0(getwd(),'/figs/image_features_reduced_3D.pdf'),width=3.5,height=2.5)

# check values are the same as in the full plot:
k_from_large <- subset(hsk, invariant == 'K' & (feature == '3D_01_target' |
                                                feature == '3D_31_target' |
                                                feature == '3D_51_target'))

summary(k_from_large$log_value)
summary(sub_dat$log_val)

# the reason for the discrepancy in the min values of log feature is that 
# the log features are computed by setting zeros to the minimum observed.
# In the full hsk data frame, that's the minimum observed over all levels.
# Other quartiles remain the same.

# "Edge density" for several features at different scales, 2D features (i.e. spatial only) --------------------- 
h <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
                            H_2D_00_target:H_2D_50_target))
h$invariant <- factor('H')
h <- melt(h, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'invariant'), variable.name = 'feature', value.name = 'value' )

# s <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
#                             S_2D_00_target:S_2D_50_target))
# s$invariant <- factor('S')
# s <- melt(s, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'invariant'), variable.name = 'feature', value.name = 'value' )

k <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
                            K_2D_00_target:K_2D_50_target))
k$invariant <- factor('K')
k <- melt(k, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'invariant'), variable.name = 'feature', value.name = 'value' )

hsk <- rbind(h, k)
hsk <- na.omit(hsk)

# remove invariant name at the front of the "feature" string so that we can compare them all at the same scale:
remove_strings <- function(i){
  string <- substring(i, 3)
  return(string)
}

# unfortunately this is super slow:
hsk$feature <- sapply(hsk$feature, remove_strings)
hsk$feature <- factor(hsk$feature)

# take log of value, correcting for zeros:
hsk$log_val <- log_convert(hsk$value)

fig <- ggplot(hsk, aes(x = log_val, y = correct, colour=invariant, fill=invariant)) + 
  facet_wrap(~ feature, ncol = 3) + 
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5)) + 
  scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
  #   geom_rug(sides = "b", alpha = 0.1) + 
  xlab('Feature value') + 
  scale_color_brewer(name='', palette='Dark2') +
  scale_fill_brewer(name='', palette='Dark2') +  
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  theme(legend.position='top')
# fig
ggsave(file=paste0(getwd(),'/figs/image_features_full_2D.pdf'),width=4.5,height=5)


# examine correlation between contrast and invariants  --------------------- 

mat <- subset(dat, sf_factor=='1.5-3', select = c(
                            K_3D_31_target, 
                            c_centre_target
                            ))
mat <- na.omit(mat)
log_mat <- mat
log_mat$K_3D_31_target <- log_convert(log_mat$K_3D_31_target)
log_mat$c_centre_target <- log_convert(log_mat$c_centre_target)

semilog_mat <- mat
semilog_mat$K_3D_31_target <- log_convert(semilog_mat$K_3D_31_target)

cor(mat)
cor(log_mat)

cor(mat, method='spearman')
cor(semilog_mat, method='spearman')

fig <- ggplot(semilog_mat, aes(x = c_centre_target, y = K_3D_31_target)) + 
  stat_binhex(bins = 50) +
#   geom_point() + 
  xlab('Contrast') + 
  ylab('Spatiotemporal feature energy') +
  scale_fill_continuous(guide=FALSE, low = "#636363", high = "#f0f0f0") +
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  theme(legend.position='top')
fig
ggsave(file=paste0(getwd(),'/figs/contrast_invariant_correlation.pdf'),width=3.5,height=3)

