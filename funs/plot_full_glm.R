## plot results of GLM fits.

rm(list=ls())
library(rstan)
library(ggplot2)
library(psybayes)
library(grid)
library(gridExtra)

source(paste0(getwd(), '/funs/helpers.R'))
source(paste0(getwd(), '/funs/plot_helpers.R'))

load(paste0(getwd(),'/output/csf_data_subset_III.RData'))
model_files <- paste0(getwd(),'/output/full_glm.RData')

load(file = model_files[1])
params <- extract(fit)

# Full glm helpers -----------------------------------
extract_params <- function(fit){
  coeffs <- c("Intercept", "Pedestal", "Increment",
              "Invariant K", "Cumulative EM",
              "0.75-1.5", "1.5-3", "3-6", "6-12", "12-24")
  
  params <- extract(fit)
  n_samples <- dim(params[[1]])[1]
  
  # extract population mean distr for each coeff:
  model_samps <- expand.grid(coeff = coeffs, sample = 1:n_samples)
  model_samps$value <- NA
  
  for (i in 1 : length(coeffs)){
    this_coeff <- coeffs[i]
    this_row <- model_samps$coeff == this_coeff
    text = paste0("this_vec <- params$mu_beta[, ", i, "]")
    eval(parse(text = text))
    # place into model samps:
    model_samps$value[this_row] <- this_vec
  }
  return(model_samps)
}

# Distributions of population mean for each coeff --------------------
samples <- extract_params(fit)

# rename intercept (because with scaled covariates, it is the intercept of the first SF):
levels(samples$coeff)[1] <- "0.375-0.75"

# plot covariates and factors separately:
slopes <- subset(samples, coeff == "Pedestal" | coeff == "Increment" | 
                   coeff == "Invariant K" | coeff == "Cumulative EM")
intercepts <- subset(samples, coeff != "Pedestal" & coeff != "Increment" & 
                       coeff != "Invariant K" & coeff != "Cumulative EM")

intercepts$coeff <- factor(intercepts$coeff) # remove extra levels

# since the intercept parameters (other than for the first level) are really *offsets* from the intercepts, to compare them on the same scale we need to add them...
for (i in 2 : length(levels(intercepts$coeff))){
  this_rows <- intercepts$coeff==levels(intercepts$coeff)[i]
  intercepts$value[this_rows] <- intercepts$value[this_rows] + 
    intercepts$value[intercepts$coeff=="0.375-0.75"]
}

fig <- ggplot(slopes, aes(x = coeff, y = value)) +
  geom_violin(size = 0.2) +
  plot_hdi_pointrange(size = 0.3) + 
  xlab("") + ylab("") +
  scale_y_continuous(name = "", breaks = c(-2, 0, 2))
fig_1 <- fig

fig <- ggplot(intercepts, aes(x = coeff, y = value)) +
  geom_violin(size = 0.2) +
  plot_hdi_pointrange(size = 0.3) + 
  xlab("") + ylab("") +
  scale_y_continuous(name = "", breaks = c(-2, 0, 2))
fig_2 <- fig


# Example model surface: Pedestal and increment --------------------

plot_dim <- 25

# generate continuous data frame for one subject, varying pedestal and increment:
# need to add some range to the predictors we don't care about so that scale() works.
surf_dat <- expand.grid(subject = 1,
                        log_ped = seq(-9, -0.5, length = plot_dim),
                        log_inc = seq(-9, -0.5, length = plot_dim),
                        sf_factor = c(levels(dat$sf_factor)[1], levels(dat$sf_factor)[2]),
                        log_K3D = seq(-2, 4, length = 3),
                        log_em_cum = seq(1, 4, length = 3))
surf_dat$sf_factor <- factor(surf_dat$sf_factor, levels = levels(dat$sf_factor))

surf_dat$p <- trial_p_glm(dat = surf_dat, param_list = params, full_model = TRUE)

# select surf dat where variables we don't want to plot are at their mean:
surf_dat <- subset(surf_dat, log_K3D == mean(log_K3D) & log_em_cum == mean(log_em_cum))

rug_correct <- subset(dat, subject == "S1" & correct == 1 & 
                        (sf_factor == levels(dat$sf_factor)[1] | sf_factor == levels(dat$sf_factor)[2]))
rug_wrong <- subset(dat, subject == "S1" & correct == 0 &
                      (sf_factor == levels(dat$sf_factor)[1] | sf_factor == levels(dat$sf_factor)[2]))

# compute percent corrects for dprimes (from ken knoblauch's dprime.mafc fun):
m <- 4
pr <- function(x, dp) dnorm(x - dp) * pnorm(x)^(m - 1)
pc_1 <- integrate(pr, lower = -Inf, upper = Inf, dp = 1)$value
pc_2 <- integrate(pr, lower = -Inf, upper = Inf, dp = 1.5)$value
pc_3 <- integrate(pr, lower = -Inf, upper = Inf, dp = 2.5)$value
dprime_breaks <- c(pc_1, pc_2, pc_3)

# do plot:
correct_col <- "#2c7bb6"
incorrect_col <- "#d7191c"

fig <- ggplot(data = surf_dat, aes(x = exp(log_ped), y = exp(log_inc), z = p)) +
  facet_wrap(~ sf_factor) +
  geom_tile(aes(fill = p)) +
  stat_contour(breaks = dprime_breaks, linetype = "dashed") +
  geom_rug(data = rug_correct, aes(x = c_centre_target, y = increment, z = NULL), 
           alpha = 0.2, colour = correct_col, sides = "bl") + 
  geom_rug(data = rug_wrong, aes(x = c_centre_target, y = increment, z = NULL), 
           alpha = 0.2, colour = incorrect_col, sides = "tr") + 
  geom_density2d(data = rug_wrong, aes(x = c_centre_target, y = increment, z = NULL), colour = incorrect_col, size = 0.3, alpha = 0.8) + 
  geom_density2d(data = rug_correct, aes(x = c_centre_target, y = increment, z = NULL), colour = correct_col, size = 0.3, alpha = 0.8) + 
  coord_cartesian(xlim = c(exp(-10), exp(0)), ylim = c(exp(-10), exp(0))) +
  scale_x_log10(name = "Pedestal contrast", breaks = c(0.001, 0.01, 0.1)) + 
  scale_y_log10(name = "Increment contrast", breaks = c(0.001, 0.01, 0.1)) + 
  scale_fill_continuous(guide = "legend", name = "Predicted P(c)", low = "#636363", high = "#f0f0f0") +
  theme(legend.position = "top")
# fig
fig_3 <- fig

# Example model surface: Pedestal and image features --------------------
surf_dat <- expand.grid(subject = 1,
                        log_ped = seq(-9, -0.5, length = plot_dim),
                        log_inc = c(-3, -2.5, -2),
                        sf_factor = c(levels(dat$sf_factor)[2]),
                        log_K3D = seq(-3, 4.5, length = plot_dim),
                        log_em_cum = c(-2, -1, 0))
surf_dat$sf_factor <- factor(surf_dat$sf_factor, levels = levels(dat$sf_factor))

surf_dat$p <- trial_p_glm(dat = surf_dat, param_list = params, full_model = TRUE)

# select surf dat where variables we don't want to plot are at their mean (i.e. scaled they will be 0):
surf_dat <- subset(surf_dat, log_em_cum == mean(log_em_cum) & log_inc == mean(log_inc))

rug_correct <- subset(dat, subject == "S1" & correct == 1 & 
                        (sf_factor == levels(dat$sf_factor)[2]))
rug_wrong <- subset(dat, subject == "S1" & correct == 0 &
                      (sf_factor == levels(dat$sf_factor)[2]))

# compute percent corrects for dprimes (from ken knoblauch's dprime.mafc fun):
m <- 4
pr <- function(x, dp) dnorm(x - dp) * pnorm(x)^(m - 1)
pc_1 <- integrate(pr, lower = -Inf, upper = Inf, dp = 1)$value
pc_2 <- integrate(pr, lower = -Inf, upper = Inf, dp = 1.5)$value
pc_3 <- integrate(pr, lower = -Inf, upper = Inf, dp = 2.5)$value
dprime_breaks <- c(pc_1, pc_2, pc_3)

# do plot:
correct_col <- "#2c7bb6"
incorrect_col <- "#d7191c"

fig <- ggplot(data = surf_dat, aes(x = exp(log_ped), y = exp(log_K3D), z = p)) +
#   facet_wrap(~ log_inc) +
  geom_tile(aes(fill = p)) + 
  stat_contour(breaks = dprime_breaks, linetype = "dashed") +
  geom_rug(data = rug_correct, aes(x = c_centre_target, y = K_3D_21_target, z = NULL), 
           alpha = 0.2, colour = correct_col, sides = "bl") + 
  geom_rug(data = rug_wrong, aes(x = c_centre_target, y = K_3D_21_target, z = NULL), 
           alpha = 0.2, colour = incorrect_col, sides = "tr") + 
  geom_density2d(data = rug_wrong, aes(x = c_centre_target, y = K_3D_21_target, z = NULL), colour = incorrect_col, size = 0.3, alpha = 0.8) + 
  geom_density2d(data = rug_correct, aes(x = c_centre_target, y = K_3D_21_target, z = NULL), colour = correct_col, size = 0.3, alpha = 0.8) + 
  scale_x_log10(name = "Pedestal contrast", breaks = c(0.001, 0.01, 0.1)) + 
  scale_y_log10(name = "Edge density", breaks = c(0.001, 0.01, 0.1, 1, 10, 100)) + 
  scale_fill_continuous(guide = "legend", name = "Predicted P(c)", low = "#636363", high = "#f0f0f0",
                        breaks = c(0.4, 0.6, 0.8)) +
  coord_cartesian(xlim = exp(c(-9.7, 0.2)), ylim = exp(c(-3.5, 5.2))) +
  theme(legend.position = "top")
# fig
fig_4 <- fig


# create a multipanel figure -----------------------------

base_size <- 8

fig_1 <- fig_1 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  ggtitle("A")

fig_2 <- fig_2 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) +
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) + 
  ggtitle("B")

fig_3 <- fig_3 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(panel.margin = unit(0.08, "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
  theme(legend.position = "top") + 
  ggtitle("C")

fig_4 <- fig_4 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(panel.margin = unit(0.08, "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
  theme(legend.position = "top") + 
  ggtitle("D")

pdf(file=paste0(getwd(),'/figs/full_glm_multipanel.pdf'),width=6,height=6)
grid.arrange(fig_1, fig_2, fig_3, fig_4, ncol = 2, widths = c(0.5, 0.5), heights = c(0.45, 0.55))
dev.off()

