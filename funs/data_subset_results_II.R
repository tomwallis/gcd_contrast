## select data subset for results part II, save to new output file.

rm(list=ls())

load(paste0(getwd(),'/output/csf_data_reduced.RData'))

# select only 1.5--3 cpd condition:
dat <- subset(dat, sf_factor == '1.5-3')
dat$sf_factor <- factor(dat$sf_factor)

# assign some crossvalidation splits:
set.seed(424242)
cv_group <- sample.int(5, size = nrow(dat), replace = TRUE)
dat$cv_group <- cv_group


# save file:
save(dat, file = paste0(getwd(),'/output/csf_data_subset.RData'))
