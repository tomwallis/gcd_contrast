# script to quickly examine various GLM forms for final model fitting.

rm(list = ls())
library(psyphy)
library(MASS)
library(psybayes)
library(ROCR)

# load data set for this section:
load(paste0(getwd(), "/output/csf_data_subset_III.RData"))
source(paste0(getwd(), "/funs/helpers.R"))

# take some logs of variables:
dat$log_K3D <- log_convert(dat$K_3D_21_target)
dat$log_ped <- log_convert(dat$c_centre_target)
dat$log_inc <- log_convert(dat$increment)
dat$log_em_cum <- log_convert(dat$em_cumDist)

# start with a big model formula. Fit to only subject 1 data:
sub_dat <- subset(dat, subject == "S1")
summary(sub_dat)

model_form <- as.formula("correct ~ scale(log_ped) + scale(log_inc) + sf_factor + stim_pos_factor + 
                         scale(log_K3D) + scale(log_em_cum)")

fit_a <- glm(data = sub_dat, formula = model_form, family = binomial())

# try stepAIC to simplify:
stepped <- stepAIC(fit_a)

new_formula <- stepped$formula
# save for future use:
# save(new_formula, file = paste0(getwd(), "/output/full_glm_formula.RData"))

# custom link:
fit <- glm(data = sub_dat, formula = new_formula, family = binomial(mafc.logit(4)), start = coef(stepped), control = list(maxit = 1000))

summary(fit, correlation = TRUE)

pairs(profile(fit))



