## create plot of TvC by binning across contrast and fitting to each level separately.

rm(list = ls())

library(ggplot2)
library(wutils)
library(psybayes)
library(psyphy)
library(plyr)

load(paste0(getwd(),'/output/csf_data_subset.RData'))

# bin pedestal contrast --------------------------

n_bins <- 6
breaks <- quantile(dat$c_centre_target, probs = seq(0, 1, length = n_bins + 1))

# where are the centres of the bins?
ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
centres <- ma(breaks)[-length(breaks)]

dat$pedestal_factor <- cut(dat$c_centre_target, breaks = breaks, include.lowest = TRUE)
dat$pedestal_centre <- NA

for (i in 1 : length(levels(dat$pedestal_factor)))
  dat$pedestal_centre[dat$pedestal_factor==levels(dat$pedestal_factor)[i]] <- centres[i]


# for each subject at each pedestal bin, fit a psychometric function as a function of increment ---------------

# return percent correct for a given dprime:
pc <- function(dp, m){
  pr <- function(x) dnorm(x - dp) * pnorm(x)^(m - 1)
  return(integrate(pr, lower=-Inf, upper = Inf)$value)
} 

t_cont <- function(pc, fit){
  eta <- fit$family$linkfun(pc)
  t_cont <- (eta - fit$coefficients[1]) / fit$coefficients[2]
  names(t_cont) <- NULL
  return(t_cont)
}

# dprime.mAFC(pc(1, 4), 4)
target_dprime = 1
threshold_pc <- pc(target_dprime, 4)

# function width (10% to 90% of function range):
width_vals <- c(0.1, 0.9)
width_vals <- 0.25 + (1 - 0.25) * width_vals

thresh_dat <- expand.grid(subject = levels(dat$subject), pedestal_factor = levels(dat$pedestal_factor))
thresh_dat$pedestal_centre <- NA
thresh_dat$y <- NA
thresh_dat$ymin <- NA
thresh_dat$ymax <- NA

for (i in 1 : length(levels(dat$pedestal_factor)))
  thresh_dat$pedestal_centre[thresh_dat$pedestal_factor==levels(thresh_dat$pedestal_factor)[i]] <- centres[i]



for (i in 1 : length(levels(dat$subject))){
  for (j in 1 : length(levels(dat$pedestal_factor))){
    subj <- levels(dat$subject)[i]
    ped <- levels(dat$pedestal_factor)[j]
    
    x <- dat$increment[dat$subject==subj & dat$pedestal_factor==ped]
    y <- dat$correct[dat$subject==subj & dat$pedestal_factor==ped]
    
    df <- data.frame(x = log(x), y = y)
    df$x[is.infinite(df$x)] <- min(df$x[!is.infinite(df$x)])
    
    # since the mafc.logit link optimiser is sensitive to starting values, initialise from normal glm parameter values:
    fit_1 <- glm(y ~ x, data = df, 
                 family = binomial())    
    
    fit <- glm(y ~ x, data = df, 
               family = binomial(mafc.logit(4)),
               start = coefficients(fit_1), maxit = 1000)
    
    # predict x for plotting:
    pred_x <- seq(log(0.00001), log(0.5), l = 200)
    pred_y <- predict(fit, newdata = data.frame(x = pred_x), type = "response")
    
    threshold <- t_cont(threshold_pc, fit)
    ymin <- t_cont(width_vals[1], fit)
    ymax <- t_cont(width_vals[2], fit)

    binned_dat <- bern_bin(df, breaks = 6, spacing = "quantile")
#     print(qplot(pred_x, pred_y, geom = "line") +  geom_pointrange(data = binned_dat, aes(x = xmid, y = ymid, ymin = ymin, ymax = ymax)) 
#           + geom_text(aes(x = -5, y = 0.9, label = round(threshold, digits = 3))) )
    
    thresh_dat$y[thresh_dat$subject==subj & thresh_dat$pedestal_factor==ped] <- exp(threshold)
    thresh_dat$ymin[thresh_dat$subject==subj & thresh_dat$pedestal_factor==ped] <- exp(ymin)
    thresh_dat$ymax[thresh_dat$subject==subj & thresh_dat$pedestal_factor==ped] <- exp(ymax)
  }
}

# plot TvC ---------------

fig <- ggplot(data = thresh_dat, aes(x = pedestal_centre, y = y, ymin = ymin, ymax = ymax))
fig <- fig + geom_pointrange()
fig <- fig + facet_wrap(~ subject, ncol = 2)
fig <- fig + scale_x_log10(name = "Pedestal contrast", breaks = c(0.01, 0.02, 0.06)) + scale_y_log10(name = "Threshold")
# fig
fig <- fig + theme_gray(base_size = 11)
ggsave(file=paste0(getwd(),'/figs/tvc_basic.pdf'),width=3.5,height=5)

