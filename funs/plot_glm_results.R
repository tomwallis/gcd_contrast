## plot results of GLM fits.

rm(list=ls())
library(rstan)
library(ggplot2)
library(psybayes)
library(grid)
library(gridExtra)

source(paste0(getwd(), '/funs/helpers.R'))
source(paste0(getwd(), '/funs/plot_helpers.R'))

load(paste0(getwd(),'/output/csf_data_subset.RData'))

model_files <- c(paste0(getwd(),'/output/single_level_glm.RData'),
                 paste0(getwd(),'/output/multilevel_glm.RData'))

# distributions for each subject's parameters -----------------------------

fig_names <- c("A", "B")

for (model in 1 : length(model_files)){
  load(file = model_files[model])
  this_model_num <- model
  samples <- param_distribution(fit = fit, dat = dat, model = this_model_num)
  
  levels(samples$type) <- c("Pr", "Post")
  levels(samples$coeff) <- c("Intercept", "Pedestal", "Increment")
  
  # this model's plot:
  fig <- ggplot(samples, aes(x = type, y = value)) +
    facet_grid(coeff ~ subject, scales = "free") +
    geom_violin(size = 0.2) +
    plot_hdi_pointrange(size = 0.3) + 
    xlab("") + ylab("") +
    scale_y_continuous(name = "", breaks = c(-2, 0, 2)) +
    ggtitle(fig_names[model]) + 
    coord_cartesian(ylim = c(-3.2, 3.2))  # zoom plot into useful range of posterior
  
  text <- paste0("fig_", model, " <- fig")
  eval(parse(text = text))
}

# example model surfaces A and B -----------------------------
surf_frame <- data.frame()
for (model in 1 : length(model_files)){
  load(file = model_files[model])
  this_model_num <- model
  this_frame <- surface_prediction(fit = fit, dat = dat, model_num = this_model_num)
  surf_frame <- rbind(surf_frame, this_frame)
}

rug_correct <- subset(dat, subject == "S1" & correct == 1)
rug_wrong <- subset(dat, subject == "S1" & correct == 0)

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

fig <- ggplot(data = surf_frame, aes(x = c_centre_target, y = increment, z = p)) +
  facet_wrap(~ model) +
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
  theme(legend.position = "top") +
  ggtitle("C")
# fig
fig_3 <- fig

# predictions for performance versus alpha plot -----------------------------

# create bins of pedestal contrast as in plot_experimental_params:
equalBreaks <- quantile(dat$c_centre_target,probs=seq(0,1,length=3))
dat$c_centre_target_factor <- cut(dat$c_centre_target,breaks=equalBreaks,include.lowest=TRUE)
levels(dat$c_centre_target_factor) <- c('Pedestal low','Pedestal high')

alpha_frame <- data.frame()
for (model in 1 : length(model_files)){
  load(file = model_files[model])
  this_model_num <- model
  this_frame <- alpha_prediction(fit = fit, dat = dat, 
                                 model_num = this_model_num,
                                 breaks = equalBreaks)
  alpha_frame <- rbind(alpha_frame, this_frame)
}

# drop incomplete cases from alpha_frame to make plots continuous:
alpha_frame <- subset(alpha_frame, complete.cases(alpha_frame))

bin_dat <- bern_bin(dat, x_name = "c_centre_target", y_name = "correct",
                    additional_factors = c("alpha", "subject"),
                    breaks = 2, spacing = "quantile", rule_of_succession = TRUE)
levels(bin_dat$c_centre_target_factor) <- c('Pedestal low','Pedestal high')

# remove cells with no trials:
bin_dat <- subset(bin_dat, n_trials > 0)

# do plot:
fig <- ggplot(bin_dat, aes(x=alpha, y=ymid)) +
  facet_grid(subject ~ c_centre_target_factor) +
  geom_ribbon(data = alpha_frame, aes(y = y, ymin=ymin, ymax=ymax, fill = model), alpha = 0.5) +
  geom_line(data = alpha_frame, aes(y = y, colour = model, fill = model), size = 0.8) +
  geom_pointrange(aes(ymin=ymin,ymax=ymax,y=ymid), colour = "black", fill = "black", size = 0.4) +
  scale_color_brewer(name = "", type = "qual", palette = 2) +
  scale_fill_brewer(name = "", type = "qual", palette = 2) +
  scale_y_continuous(name='Proportion correct',limits=c(0,1),breaks=seq(0,1,l=3)) + 
  scale_x_continuous(name='Multiplication factor (alpha)', breaks=seq(2, 6, l = 3))
# fig
fig_4 <- fig + ggtitle("D")

# create a multipanel figure -----------------------------

base_size <- 8

fig_1 <- fig_1 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))

fig_2 <- fig_2 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) +
  theme(strip.background = element_rect(fill = "grey80", colour = "White"))

fig_3 <- fig_3 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(panel.margin = unit(0.08, "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
  theme(legend.position = "top")
#   theme(legend.margin = unit(0.01, "cm"), legend.key.size = unit(0.02, "cm") )

fig_4 <- fig_4 + theme_minimal(base_size=base_size) + 
  theme(plot.margin =unit(rep(0.02, times = 4), "inches") ) + 
  theme(panel.margin = unit(0.08, "inches") ) + 
  theme(strip.background = element_rect(fill = "grey80", colour = "White")) +
  theme(legend.position = "top")

pdf(file=paste0(getwd(),'/figs/glm_multipanel.pdf'),width=6.5,height=6)
grid.arrange(fig_1, fig_2, fig_3, fig_4, ncol = 2, widths = c(0.5, 0.5), heights = c(0.45, 0.55))
dev.off()

