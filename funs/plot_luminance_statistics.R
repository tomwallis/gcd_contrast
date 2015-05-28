## plot contrast stats:

# plot distribution of contrast observed in contrast increment detection experiments.
# TSAW.

rm(list = ls())

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)
library(grid)
library(gridExtra)
library(mgcv)

# load data frame:
load(paste(getwd(),'/output/csf_data_full.RData',sep=''))

# # plot histogram:
# fig <- ggplot(dat,aes(x=lum_target)) + 
# #   facet_wrap(~ sf_factor, ncol=2, scales='free_y') + 
#   geom_histogram(aes(colour=NULL),position='identity') + 
#   scale_y_continuous(name='',breaks=NULL) + 
#   xlab('Pixel intensity')
# # ggsave(file=paste(getwd(),'/figs/contrast_histogram.pdf',sep=""),width=3.5,height=3)
# 
# fig_1 <- fig + ggtitle("A") + 
#   theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
#   theme_minimal(base_size=8) + 
#   theme(strip.background = element_rect(fill = "grey80", colour = "White"))

# how many times do we observe zero luminance?
zeros <- dat$c_centre_target[dat$lum_target==0]
print(n_zeros <- length(zeros))


# plot correlation with contrast, 1.5--3 cpd -----------------------------


# mat <- subset(dat, sf_factor=='1.5-3', select = c(c_centre_target,
#                                                   lum_target))

mat <- subset(dat, select = c(c_centre_target,
                              lum_target))

mat <- na.omit(mat)
cor(mat)
cor(mat, method='spearman')

fig <- ggplot(mat, aes(x = lum_target, y = c_centre_target)) + 
  stat_binhex(bins = 40) +
#   stat_density2d() +
  #   geom_point() + 
  xlab('Mean pixel intensity') + 
  ylab('Contrast') +
  scale_fill_continuous(guide=FALSE, low = "#636363", high = "#f0f0f0") +
  theme(legend.position='top')
# fig

fig_2 <- fig + ggtitle("A") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))


# plot relationship with performance ----

# for three bins of contrast
dat$contrast_bin <- cut(dat$c_centre_target, breaks=c(0, 0.01, 0.05, 0.2))

print(levels(dat$contrast_bin)[2])
sub_dat <- subset(dat, contrast_bin==levels(dat$contrast_bin)[2])

# correlation in our subset:
cor(sub_dat$lum_target, sub_dat$c_centre_target, method='spearman', use='complete.obs')


fig <- ggplot(dat, aes(x=lum_target, y=correct)) +
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5), colour='black') +
  stat_smooth(data=sub_dat, 
              family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5), 
              colour='black', linetype='dashed') +
  geom_rug(sides = "b", alpha = 0.1) + 
  scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
  xlab('Mean pixel intensity') + 
  theme(legend.position='top')

fig_3 <- fig + ggtitle("B") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))  
  


# check that within the subset, contrast and performance are related:
fig <- ggplot(sub_dat, aes(x=c_centre_target, y=correct)) +
  stat_smooth(data=sub_dat, 
              family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5), 
              colour='black', linetype='dashed') +
  geom_rug(sides = "b", alpha = 0.1) + 
  scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
  xlab('Contrast') + 
  theme(legend.position='top')

fig_4 <- fig + ggtitle("C") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))  

# create a multipanel figure:
pdf(file=paste0(getwd(),'/figs/luminance_statistics.pdf'),width=6,height=2.5)

# grid.arrange(fig_1, fig_2, fig_3, ncol = 2, widths = c(1, 1))
grid.arrange(fig_2, fig_3, fig_4, ncol = 3, widths = c(1, 1))

dev.off()


# # examine relationship between luminance and performance  --------------------- 
# ## Idea for plot: could do panels of S, W, E and N, then plot separate lines for target positions -- shows long-range relationships.
# 
# sub_dat <- subset(dat, select = c(unique_id,subject,sf_factor,correct,
#                                   lum_S:lum_N, stim_pos_factor))
# sub_dat <- na.omit(sub_dat)
# 
# sub_dat <- melt(sub_dat, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct', 'stim_pos_factor'), variable.name = 'feature', value.name = 'value' )
# levels(sub_dat$feature)
# 
# fig <- ggplot(sub_dat, aes(x = value, y = correct, colour = stim_pos_factor)) + 
#   facet_wrap(~ feature, ncol = 2) + 
#   stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5),aes(fill=NULL)) + 
#   geom_rug(sides = "b", alpha = 0.1) + 
#   scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
#   xlab('Luminance') + 
#   theme_minimal(base_size=11) + 
#   theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
#   theme(legend.position='top')
# # fig
# ggsave(file=paste0(getwd(),'/figs/image_features_luminance.pdf'),width=4,height=4)



