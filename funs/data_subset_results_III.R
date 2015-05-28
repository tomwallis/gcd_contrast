## select data subset for results part II, save to new output file.

rm(list=ls())

load(paste0(getwd(),'/output/csf_data_full.RData'))
source(paste0(getwd(),"/funs/helpers.R"))

# exclude TW, AJ and LB (few trials)
dat <- dat[dat$subject!='AJ',]
dat <- dat[dat$subject!='LB',]
dat <- dat[dat$subject!='TSAW',]
dat$subject <- factor(dat$subject) # get rid of empty levels.

# rename subjects to anonymous codes:
levels(dat$subject) <- c('S1','S2','S3','S4','S5') # corresponds to MACD, PJB, LAL, CPT and AM.


# select variables to include in GLM (based on exploratory plots):
dat <- subset(dat, select = c(unique_id : resp_factor, 
                              c_centre_target, c_surround_target,
                              increment : alpha_factor,
                              em_nextSacc_dir_relative,
                              em_cumDist,
                              K_3D_21_target))


# drop missing cases (ouch!)
dat <- na.omit(dat)
dat$subject <- factor(dat$subject) # get rid of empty levels.


# assign test / training split:
set.seed(424242)
dat$cv_group <- "train"
test_rows <- sample.int(nrow(dat), size = nrow(dat) * 0.2)
dat$cv_group[test_rows] <- "test"
dat$cv_group <- factor(dat$cv_group)

# check that each subject and sf is represented in the test set:
library(dplyr)
n_trials <- dat %.%
  group_by(subject, sf_factor, cv_group) %.%
  summarise(n_trials = length(correct))
n_trials

# save file:
save(dat, file = paste0(getwd(),'/output/csf_data_subset_III.RData'))
