## Master file for producing analysis. Calls functions in /funs, prints figures to /figs.

# Import and clean data files -------------------------------------------
rm(list = ls())
source(paste0(getwd(),'/funs/import_data.R'))
# create a full data set with all predictors, subjects and missing values, 
# as well as a reduced set for contrast analysis.
import_data()

# examine which subjects did which alpha levels -------------------------
load(paste0(getwd(),'/output/csf_data_reduced.RData'))
library(plyr)
sum_dat <- ddply(dat,.(alpha,subject),summarise,count=length(correct))


# Results part I ---------------------------
# Create data summary table ---------------------------------------------
source(paste0(getwd(),'/funs/summary_table.R'))

# Do primary data plots --------------------------------------------------
source(paste0(getwd(),'/funs/plot_contrast_statistics.R'))
source(paste0(getwd(),'/funs/plot_luminance_statistics.R'))
source(paste0(getwd(),'/funs/plot_experimental_params.R'))
source(paste0(getwd(),'/funs/plot_eye_movements.R'))
source(paste0(getwd(),'/funs/plot_image_features.R'))


# Results part II ---------------------------
source(paste0(getwd(),'/funs/data_subset_results_II.R'))
source(paste0(getwd(),'/funs/tvc_plot_psyphy.R'))

# Sampling operations (take time) ---------------------------------
# it might be an idea to run these on a cluster, or at least wait a while:
source(paste0(getwd(),'/funs/fit_single_level_transducer.R'))
source(paste0(getwd(),'/funs/fit_single_level_transducer_B.R'))
source(paste0(getwd(),'/funs/fit_single_level_glm.R'))
source(paste0(getwd(),'/funs/fit_multilevel_glm.R'))

source(paste0(getwd(),'/funs/cv_model_1.R'))
source(paste0(getwd(),'/funs/cv_model_2.R'))
source(paste0(getwd(),'/funs/cv_model_3.R'))
source(paste0(getwd(),'/funs/cv_model_4.R'))

# Plotting results from sampling ---------------------------------------

source(paste0(getwd(),'/funs/produce_mcmc_diagnostics.R'))

# plot results from model fitting:
source(paste0(getwd(),'/funs/plot_glm_results.R'))
source(paste0(getwd(),'/funs/plot_transducer_results.R'))
source(paste0(getwd(),'/funs/plot_correlations.R'))

# plot crossvalidation results:
source(paste0(getwd(),'/funs/cv_mean_performance.R'))
source(paste0(getwd(),'/funs/cv_max_likelihood.R'))
source(paste0(getwd(),'/funs/plot_crossvalidation_results.R'))

# Results part III ---------------------------
source(paste0(getwd(),'/funs/data_subset_results_III.R'))
source(paste0(getwd(),'/funs/glm_explore.R'))

source(paste0(getwd(),'/funs/fit_full_glm.R')) # takes a long time...

source(paste0(getwd(),'/funs/mcmc_diagnostics_full_glm.R'))
source(paste0(getwd(),'/funs/plot_full_glm.R'))
source(paste0(getwd(),'/funs/do_full_glm_crossval.R'))

