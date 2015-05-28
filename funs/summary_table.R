## Print a table of trials by subject.

rm(list = ls())

library(plyr)
library(xtable)
load(paste(getwd(),'/output/csf_data_reduced.RData',sep=''))

sum_dat <- ddply(dat,.(subject,sf_factor),summarise,N=length(correct))

table_dat <- subset(dat,select=c("subject","sf_factor"))

tab <- table(table_dat)
tab <- addmargins(tab,margin=c(1,2))

print_tab <- xtable(tab,digits=0,
                     label='tab:trial_numbers_nas',
                     caption='Number of trials for all subjects at each target spatial frequency band.')

write(print(print_tab), file = paste0(getwd(),'/figs/summary_table.tex'))
