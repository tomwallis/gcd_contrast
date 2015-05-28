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

# load data frame:
load(paste(getwd(),'/output/csf_data_reduced.RData',sep=''))
source(paste0(getwd(), '/funs/helpers.R'))

# plot histogram:
fig <- ggplot(dat,aes(x=c_centre_target)) + 
  facet_wrap(~ sf_factor, ncol=2, scales='free_y') + 
  geom_histogram(aes(colour=NULL),position='identity',binwidth=.005) + 
  scale_x_continuous(breaks = c(0, 0.1, 0.2)) + 
  scale_y_continuous(name='',breaks=NULL) + 
  xlab('Contrast Energy')
# ggsave(file=paste(getwd(),'/figs/contrast_histogram.pdf',sep=""),width=3.5,height=3)

fig_1 <- fig + ggtitle("A") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))

# how many times do we observe zero contrast?
zeros <- dat$c_centre_target[dat$c_centre_target==0]
print(n_zeros <- length(zeros))

# show a log-scaled contrast histogram too:
dat$log_contrast = log_convert(dat$c_centre_target)

fig <- ggplot(dat,aes(x=log_contrast)) + 
  facet_wrap(~ sf_factor, ncol=2, scales='free_y') + 
  geom_histogram(aes(colour=NULL),position='identity', binwidth=0.2) + 
#   scale_x_continuous(breaks = c(0, 0.1, 0.2)) + 
  scale_y_continuous(name='',breaks=NULL) + 
  xlab('Contrast Energy [log units]')
ggsave(file=paste0(getwd(),'/figs/log_contrast_histogram.png'))

# plot contrast correlations -----------------------------

sub_dat <- subset(dat,select=c(c_centre_band_0_target:c_centre_band_5_target))
# rename columns:
colnames(sub_dat) <- as.character(levels(dat$sf_factor))
cor_frame <- melt(cor(sub_dat, method = "spearman"))
cor_frame$Var1 <- factor(cor_frame$Var1,levels(dat$sf_factor))
cor_frame$Var2 <- factor(cor_frame$Var2,levels(dat$sf_factor))


fig <- ggplot(data = cor_frame, aes(x=Var1, y=Var2, fill=value, label=round(value, 2))) + 
  scale_fill_gradient() + geom_tile() + geom_text(size=2.5) + 
  xlab('') + ylab('') + 
  scale_fill_gradient2(name='Correlation',limits=c(-1,1), guide=FALSE,
                       low = muted("red"),
                       mid = "white", high = muted("blue"))
# ggsave(fig,file=paste0(getwd(),'/figs/contrast_correlations.pdf'),width=3.5,height=3.5)

fig_2 <- fig + ggtitle("B") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
  theme(axis.line = element_blank(), axis.text.x = element_text(angle=45)) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))

# plot overlap of spatial frequency bands (reviewer request) ----

# the data file was sent to me by Michael, email dated 2015.02.08
band_output <- readRDS(paste0(getwd(), '/output/band_output.rds'))
# band_output <- load(paste0(getwd(), '/output/band_output.RData'))

# rearrange for plotting:
bands <- melt(band_output, id.vars='x', variable.name='band', value.name='response')
levels(bands$band) <- rev(levels(dat$sf_factor))

# reverse levels for plotting:
bands$band <- factor(bands$band, levels=rev(levels(bands$band)))


# colorbrewer dark2 palette:
pal <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02')


# plot:
fig <- ggplot(bands, aes(x=x, y=response, colour=band)) +
  geom_line() + 
  scale_x_log10(name='Spatial frequency (cpd)') + 
  scale_y_continuous(name='Normalised filter response', limits=c(0, 1), breaks=seq(0, 1, length=3)) + 
  scale_color_manual(name='', values=pal)
# fig

fig_3 <- fig + ggtitle("C") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
#   theme(axis.line = element_blank(), axis.text.x = element_text(angle=45)) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  theme(legend.position='top')

  
# plot falloff in correlation for some spatial scales (reviewer request) ----

# need to express some distance from each filter to every other, then access the 
# correlations in cor_frame. Will hack this for now.

delta_correlations <- data.frame()
n_filts <- 6
filt_nums <- 1:n_filts
for (i in 1:n_filts){
  i_filt <- levels(dat$sf_factor)[i]
  delta <- filt_nums - i
  for (j in 1:n_filts){
    j_filt <- levels(dat$sf_factor)[j]
    # access correlation coefficient i,j
    corr <- cor_frame$value[cor_frame$Var1==i_filt & cor_frame$Var2==j_filt]
    
    this_df <- data.frame(delta=delta[j], filter=i_filt, correlation=corr)
    delta_correlations <- rbind(delta_correlations, this_df)
  }
}

# throw out some levels to declutter the plot:
delta_correlations <- subset(delta_correlations, filter=='.375-.75' | 
                               filter=='1.5-3' |
                               filter=='12-24')

fig_4 <- ggplot(delta_correlations, aes(x=delta, y=correlation, colour=filter)) +
  geom_point() + 
  geom_line() + 
  scale_x_continuous(name='Filter level difference (lower = more coarse)', breaks=seq(-5, 5, l=11)) + 
  scale_y_continuous(name='Correlation', limits=c(0, 1), breaks=seq(0, 1, length=3)) + 
  scale_color_manual(name='', values=c(pal[1], pal[3], pal[6]))

fig_4 <- fig_4 + ggtitle("D") + 
  theme(plot.margin =unit(rep(0.05, times = 4), "inches") ) + 
  theme_minimal(base_size=8) + 
#   theme(axis.line = element_blank(), axis.text.x = element_text(angle=45)) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  theme(legend.position='top')


# fig_4            

# blank panel (from http://www.r-bloggers.com/extra-extra-get-your-gridextra/):

blank_panel <- grid.rect(gp=gpar(col="white"))
# df <- data.frame()
# blank_panel <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100) + ggtitle("C") + theme.blank()
# blank_panel

# create a multipanel figure:
pdf(file=paste0(getwd(),'/figs/contrast_statistics.pdf'),width=5.5,height=5.5)

grid.arrange(fig_1, fig_2, fig_3, fig_4, ncol = 2, widths = c(1.2, 1))

dev.off()

