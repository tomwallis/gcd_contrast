## Explore data via ggplot smooths.
# TSAW

rm(list = ls())

library(ggplot2)
library(mgcv)
library(plyr)
library(wutils)
library(psybayes)
library(reshape2)
library(grid)

load(paste0(getwd(),'/output/csf_data_reduced.RData'))
source(paste0(getwd(), '/funs/helpers.R'))


# Performance as a function of alpha, pedestal --------------------- 

# create bins of pedestal contrast:
equalBreaks <- quantile(dat$c_centre_target,probs=seq(0,1,length=3))
dat$c_centre_target_factor <- cut(dat$c_centre_target,breaks=equalBreaks,include.lowest=TRUE)
levels(dat$c_centre_target_factor) <- c('low','high')

bin_dat <- bern_bin(dat, x_name = "c_centre_target", y_name = "correct",
                 additional_factors = c("alpha", "subject", "sf_factor"),
                 breaks = 2, spacing = "quantile", rule_of_succession = TRUE)
levels(bin_dat$c_centre_target_factor) <- c('low','high')

# remove cells with no trials:
bin_dat <- subset(bin_dat, n_trials > 0)

fig <- ggplot(dat,aes(x=alpha,y=correct,colour=c_centre_target_factor,shape=c_centre_target_factor,fill=c_centre_target_factor)) +
  facet_grid(subject ~ sf_factor) +
  geom_hline(aes(yintercept=0.25), size=.5, linetype="dashed", colour="grey") + 
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=3),aes(fill=NULL)) + 
  geom_linerange(data=bin_dat,aes(ymin=ymin,ymax=ymax,y=ymid)) + geom_point(data=bin_dat,aes(y=ymid))
fig <- greyfig(fig,name='Pedestal contrast')
fig <- fig + scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) +
  scale_x_continuous(name='Multiplication factor (alpha)', breaks=seq(2, 6, l = 3)) + 
  theme_minimal(base_size=11) + 
  theme(legend.position='top') + 
  theme(panel.margin = unit(0.2, "inches")) +
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))
# fig
ggsave(file=paste0(getwd(),'/figs/explore_performance.pdf'),width=7,height=6.5)



# create averages across subjects ------------------------
average_dat <- ddply(bin_dat, .(c_centre_target_factor, sf_factor, alpha), summarise,
                     y = mean(ymid, na.rm = TRUE), 
                     ymax = mean(ymid, na.rm = TRUE) + sd(ymid, na.rm = TRUE), 
                     ymin = mean(ymid, na.rm = TRUE) - sd(ymid, na.rm = TRUE))
fig2 <- ggplot(average_dat,aes(x=alpha,y=correct,colour=c_centre_target_factor,shape=c_centre_target_factor,fill=c_centre_target_factor))
fig2 <- fig2 + facet_wrap(~ sf_factor, ncol = 3)
fig2 <- fig2 + geom_linerange(data=average_dat,aes(ymin=ymin,ymax=ymax,y=y)) + geom_point(data=average_dat,aes(y=y)) 
fig2 <- greyfig(fig2,name='Pedestal contrast')
fig2 <- fig2 + scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3))
fig2 <- fig2 + scale_x_continuous(name='Multiplication factor (alpha)', breaks=seq(2, 6, l = 3))
fig2 <- fig2 + theme_grey(base_size=11)
fig2 <- fig2 + theme(legend.position='top')
library(grid)
fig2 <- fig2 + theme(panel.margin = unit(0.1, "inches"))
# fig2
# ggsave(file=paste0(getwd(),'/figs/explore_performance_average.pdf'),width=3.5,height=3.5)



## Reviewer: why does performance in high pedestal case (in explore_performance plot) appear higher than low pedestal? How does this fit with masking effects in model?
## Difference is that explore_performance shows multiplication factor, not increment or pedestal.
## Let's create some plots that make this more intuitive. Focus just on S1 data.

sub_dat <- subset(dat, subject=='S1' & sf_factor=='1.5-3')
# sub_dat$log_increment <- log_convert(sub_dat$increment)
sub_dat$alpha_factor <- factor(sub_dat$alpha_factor)
low_rug <- subset(sub_dat, c_centre_target_factor=='low')
high_rug <- subset(sub_dat, c_centre_target_factor=='high')



fig <- ggplot(sub_dat,aes(x=increment,y=correct, colour=c_centre_target_factor)) +
  geom_hline(aes(yintercept=0.25), size=.5, linetype="dashed", colour="grey") + 
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5),
              aes(fill=NULL)) + 
  geom_rug(data=low_rug, sides='b') + 
  geom_rug(data=high_rug, sides='t')
fig <- fig + scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) +
  scale_x_continuous(name='Increment contrast') + 
  theme_minimal(base_size=8) +
  scale_color_manual(name = "Pedestal contrast", values=c('#969696', '#252525')) + 
  theme(legend.position='top')
fig
ggsave(file=paste0(getwd(),'/figs/increment_vs_performance.pdf'),width=3.5,height=3.5)

