# do plots of eye movement-related stuff.

rm(list = ls())

library(ggplot2)
library(mgcv)
library(plyr)
library(wutils)
library(psybayes)
library(reshape2)
library(gridExtra)
library(circular)

# needs "full" data file.
load(paste0(getwd(),'/output/csf_data_full.RData'))

# drop extra subjects:
# exclude TW, AJ and LB (few trials)
dat <- dat[dat$subject!='AJ',]
dat <- dat[dat$subject!='LB',]
dat <- dat[dat$subject!='TSAW',]
dat$subject <- factor(dat$subject) # get rid of empty levels.

# rename subjects to anonymous codes:
levels(dat$subject) <- c('S1','S2','S3','S4','S5') # corresponds to MACD, PJB, LAL, CPT and AM.


# Histograms of eye movement directions before and after target -----------------------
dir_frame <- subset(dat, select = c(unique_id, subject, sf_factor, correct, em_prevSacc_dir_relative, em_nextSacc_dir_relative))
dir_frame <- melt(dir_frame, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct'), variable.name = 'saccade', value.name = 'direction')
levels(dir_frame$saccade) <- c('Previous Saccade','Next Saccade')
# remove NAs for plotting:
dir_frame <- na.omit(dir_frame)

# move from -pi +pi range to 0 -- 2pi:
dir_frame$direction <- dir_frame$direction + pi

# estimate circular density using density.circular from the circular package:
bandwidth <- 0.02

dir <- subset(dir_frame, saccade == "Previous Saccade", select = direction)
dir <- circular(dir$direction, type = "angles", units = "radians")
circ_density <- density.circular(dir, bw = bandwidth, kernel = "wrappednormal", n = 1000)
circ_frame <- data.frame(x = circ_density$x, y = circ_density$y, saccade = "Previous Saccade")

dir <- subset(dir_frame, saccade == "Next Saccade", select = direction)
dir <- circular(dir$direction, type = "angles", units = "radians")
circ_density <- density.circular(dir, bw = bandwidth, kernel = "wrappednormal", n = 1000)
tmp_frame <- data.frame(x = circ_density$x, y = circ_density$y, saccade = "Next Saccade")

circ_frame <- rbind(circ_frame, tmp_frame)

fig <- ggplot(circ_frame, aes(x = x, y = y)) +
  facet_wrap(~ saccade, ncol = 2) +
  geom_line() +
  scale_x_continuous(name="Saccade direction relative to target (north)", 
                     breaks = c(0, pi/2, pi, 3*pi/2), labels = NULL) +
  scale_y_continuous(name="Density", breaks = c(0.1, 0.2, 0.3, 0.4)) +
  coord_polar(start=pi) +
  theme_grey(base_size=11)
# fig
# ggsave(file=paste0(getwd(),'/figs/em_direction_histograms.pdf'),width=3.5,height=3)

fig_1 <- fig + ggtitle("A")

# Histograms of eye movement amplitudes before and after target -----------------------
amp_frame <- subset(dat, select = c(unique_id, subject, sf_factor, correct, em_prevSacc_amp, em_nextSacc_amp))
amp_frame <- melt(amp_frame, id.vars = c('unique_id', 'subject', 'sf_factor', 'correct'), variable.name = 'saccade', value.name = 'amplitude')
levels(amp_frame$saccade) <- c('Previous Saccade','Next Saccade')
# remove NAs for plotting:
amp_frame <- na.omit(amp_frame)

labels <- ddply(amp_frame, .(saccade), summarise, label = paste0("Mean: ", round(mean(amplitude), digits = 2)))
labels$x <- 0.5
labels$y <- c(200, 400)

fig <- ggplot(amp_frame, aes(x = saccade, y = amplitude)) +
  geom_violin() +
  plot_hdi_pointrange() +
  xlab("") + 
  scale_y_log10(name="Eye movement amplitude (deg)", breaks = c(0.3, 1, 3, 10, 30, 50)) +
  theme_grey(base_size=11)
# fig
# ggsave(file=paste0(getwd(),'/figs/em_amplitude_histograms.pdf'),width=3.5,height=3)

fig_2 <- fig + ggtitle("B")

# Cumulative eye movement distance and proportion correct -----------------------
fig <- ggplot(dat, aes(x = em_cumDist, y = correct))
fig <- fig + stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=5), colour = "black")
fig <- fig + geom_rug(sides = "b", alpha = 0.5)
fig <- fig + scale_x_log10(name="Cumulative eye movement distance (deg)", breaks=c(3,10,30)) 
fig <- fig + scale_y_continuous(name="Proportion correct",breaks=seq(0,1,by=.5), limits=c(0,1))
fig <- fig + theme_grey(base_size=11)
#fig
# ggsave(file=paste0(getwd(),'/figs/em_cum_dist.pdf'),width=3.5,height=3)

fig_3 <- fig + ggtitle("C")

# split into factors by saccade amplitude, timing, and direction relative to target -------------------
# breaks motivated by experiment spatial coord (2 deg eccent, sd = 1 deg):
spatial_breaks <- c(0,1,3,Inf)
temporal_breaks <- c(0,120,480,Inf)

dat$amplitude_factor <- cut(dat$em_nextSacc_amp,breaks=spatial_breaks,include.lowest=TRUE)
levels(dat$amplitude_factor)[length(levels(dat$amplitude_factor))] <- '(3,max]'

dat$time_factor <- cut(dat$em_nextSacc_time_StimOnToOn,breaks=temporal_breaks,include.lowest=TRUE)
levels(dat$time_factor)[length(levels(dat$time_factor))] <- '(480,max]'

# drop NAs for plotting:
dat <- subset(dat,time_factor != "NA's")

# plot as polar plot:

fig <- ggplot(dat,aes(x=em_nextSacc_dir_relative,y=correct)) 
knots <- seq(-pi,pi, by=pi/4)
fig <- fig + stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cc",k=length(knots)),knots=list(x=knots),colour="black",size=0.5,fill="grey50")
fig <- fig + facet_grid( amplitude_factor ~ time_factor)
fig <- fig + scale_x_continuous(name="Eye movement direction relative to target",
                                breaks=c(-pi, -pi/2, 0, pi/2, pi), labels = c('A','L','T','R','A')) 
fig <- fig + scale_y_continuous(name="Proportion correct",breaks=seq(0,1,by=.5), limits=c(0,1))
# fig <- fig + coord_polar(start=pi)
fig <- fig + theme_grey(base_size=11)
library(grid)
fig <- fig + theme(panel.margin = unit(0.15, "inches"))
# fig
# ggsave(file=paste0(getwd(),'/figs/em_performance.pdf'),width=7,height=6)

fig_4 <- fig + ggtitle("D")

# don't facet by amplitude and time?
fig <- ggplot(dat,aes(x=em_nextSacc_dir_relative,y=correct)) 
knots <- seq(-pi,pi, by=pi/4)
fig <- fig + stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cc",k=length(knots)),knots=list(x=knots),colour="black",size=0.5,fill="grey80")
fig <- fig + scale_x_continuous(name="Eye movement direction relative to target",breaks=c(-pi, -pi/2, 0, pi/2, pi), labels = c('Away','Left','Towards','Right', 'Away')) 
fig <- fig + scale_y_continuous(name="Proportion correct",breaks=seq(0,1,by=.5), limits=c(0,1))
# fig <- fig + coord_polar(start=pi)
fig <- fig + theme_grey(base_size=11)
# fig
# ggsave(file=paste0(getwd(),'/figs/em_performance_not_split.pdf'),width=3.5,height=3)


# Multipanel eye movements figure -----------------
base_size <- 8

fig_1 <- fig_1 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))

fig_2 <- fig_2 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) +
  theme(legend.position = "top") + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))
  
fig_3 <- fig_3 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))

fig_4 <- fig_4 + theme_minimal(base_size=base_size) +
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(panel.margin = unit(0.05, "inches")) +
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))


pdf(file=paste0(getwd(),'/figs/em_multipanel.pdf'),width=6.5,height=5.5)
grid.arrange(fig_1, fig_2, fig_3, fig_4, ncol = 2, widths = c(0.5, 0.5))
dev.off()


# Sensitivity as a function of timing relative to previous / next saccade ----------
prev_timing <- subset(dat, select = c(unique_id, subject, sf_factor, correct, em_prevSacc_time_StimOnToOff, em_prevSacc_dir_relative))
next_timing <- subset(dat, select = c(unique_id, subject, sf_factor, correct, em_nextSacc_time_StimOnToOn, em_nextSacc_dir_relative))

time_cutoff <- 2000

prev_timing <- subset(prev_timing, em_prevSacc_time_StimOnToOff > -time_cutoff)
next_timing <- subset(next_timing, em_nextSacc_time_StimOnToOn < time_cutoff)

# create factors from "towards target" and "other":
rad_cut <- 22.5 * (pi / 180)

prev_timing$dir_factor <- "Other"
prev_timing$dir_factor[prev_timing$em_prevSacc_dir_relative >= -rad_cut & 
                         prev_timing$em_prevSacc_dir_relative <= rad_cut] <- "Towards"

next_timing$dir_factor <- "Other"
next_timing$dir_factor[next_timing$em_nextSacc_dir_relative >= -rad_cut & 
                         next_timing$em_nextSacc_dir_relative <= rad_cut] <- "Towards"


prev_timing$dir_factor <- factor(prev_timing$dir_factor)
next_timing$dir_factor <- factor(next_timing$dir_factor)

fig <- ggplot(prev_timing, aes(x = em_prevSacc_time_StimOnToOff, y = correct, linetype = dir_factor)) + 
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=8), colour = "black") + 
#   geom_rug(sides = "b", alpha = 0.5) + 
  scale_x_continuous(name="Time from offset of previous saccade (ms)") +
  scale_y_continuous(name="Proportion correct",breaks=seq(0,1,by=.5), limits=c(0,1)) + 
  scale_linetype_discrete(name="") +
  theme_minimal(base_size=base_size)
fig_5 <- fig + ggtitle("A")


fig <- ggplot(next_timing, aes(x = em_nextSacc_time_StimOnToOn, y = correct, linetype = dir_factor)) +
#   facet_grid(sf_factor ~ subject) +
  stat_smooth(family=binomial(),method="gam",formula=y ~ s(x,bs="cs",k=8), colour = "black") + 
#   geom_rug(sides = "b", alpha = 0.5) + 
  scale_x_continuous(name="Time to onset of next saccade (ms)") +
  scale_y_continuous(name="Proportion correct",breaks=seq(0,1,by=.5), limits=c(0,1)) + 
  scale_linetype_discrete(name="") +
  theme_minimal(base_size=base_size) 
# fig
fig_6 <- fig + ggtitle("B")

pdf(file=paste0(getwd(),'/figs/em_timing.pdf'),width=3.5,height=4)
grid.arrange(fig_5, fig_6, ncol = 1, widths = c(0.5, 0.5))
dev.off()


# Geotopic responses plot -------------------------
# Plot responses relative to target location for large saccades made during target presentation as in Dorr & Bex...
# select saccades made towards the target (central +- 45 degrees i.e. 0.785 radians):

# I think this will require a new timing variable from Michael (time from stim off to saccade off). I don't think
# this can be constructed out of any of the current timings.

# just for now ignore timing, to create correct groups.
# geoFrame <- subset(dat,em_nextSacc_dir_relative>=-0.785 & em_nextSacc_dir_relative<=0.785)
# 
# geoFrame$congruency <- 'orthogonal'
# geoFrame$congruency <- factor(geoFrame$congruency,levels=c('orthogonal','correct','incongruent'))
# # arrange choices similar to Dorr and Bex (note they used eye movement direction relative to target):
# for (i in 1:length(unique(geoFrame$stim_pos))){
#   thisPos <- unique(geoFrame$stim_pos)[i]
#   if(thisPos==1) antiPos <- 9
#   if(thisPos==3) antiPos <- 7
#   if(thisPos==7) antiPos <- 3
#   if(thisPos==9) antiPos <- 1
#   geoFrame$congruency[geoFrame$stim_pos==thisPos & geoFrame$resp==thisPos] <- 'correct'
#   geoFrame$congruency[geoFrame$stim_pos==thisPos & geoFrame$resp==antiPos] <- 'incongruent'
# }
# summary(geoFrame)