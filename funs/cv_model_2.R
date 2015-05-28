# script to do crossvalidation for one model. Separate script for easier calling on cluster.

rm(list=ls())
source(paste0(getwd(),'/funs/crossval_functions.R'))
do_crossval(model = 2)
