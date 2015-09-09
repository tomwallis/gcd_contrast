# Detection of gaze-contingent contrast increments in natural image movies

These are R scripts to reproduce the analyses in Wallis, Dorr & Bex (in press).

These files are released under the GPL-3 License. I don't offer any support. You are welcome to send me an email but I make no claim that I will (be able to) help, particularly if you are attempting to use a different platform or setup than mine (I'm running on OSX 10.9 with RStudio; the MCMC sampling was conducted on a Debian distribution).

**IMPORTANT**: If you adapt any part of these scripts for your own academic work, please help keep me employed by citing the following paper:

Wallis, T. S. A., Dorr, M., & Bex, P. J. (2015). Sensitivity to gaze-contingent contrast increments in naturalistic movies: An exploratory report and model comparison. *Journal of Vision, 15*(8), 3. http://doi.org/10.1167/15.8.3

## Subdirectories in this repository

The main file is `analysis_master.R` in the parent directory. Use this as a guide to follow all the analyses.

Other directories are:
  * /output/ contains data used for analysis. These are currently distributed as R binaries to reduce file size, but if you want .csv files let me know.
  * /funs/ contains all scripts and functions to produce the analysis and figures. There are also a bunch of legacy functions included for completeness that were not used in the final paper, and may no longer run.


## Very raw data

We have chosen not to provide the raw-est of raw data or MCMC samples here due to size limitations (these files are very large). If you need them, please contact me. They may go up on a different venue at some time; if so I will provide a link here.


## Reproducing the analysis

The analyses in the paper were performed in R using various libraries. By looking through the `analysis_master.R` file you can determine which scripts in the /funs/ directory to run in which order (or just source `analysis_master.R`. There may be problems with paths, particularly if you're running a windows system (this code uses hard-coded `/` characters for directory separation; sorry).

## Dependencies

I've tried to compile a list of dependencies for these scripts here. These packages are required to run all the scripts to reproduce the analysis. If I've left something off please let me know.

  * Rstan (available from www.mc-stan.org)
  * wutils (available from my github [here](https://github.com/tomwallis/wutils))
  * psybayes (available from my github [here](https://github.com/tomwallis/psybayes))

Everything else is on CRAN:

  * plyr
  * dplyr
  * xtable
  * ggplot2
  * reshape2
  * scales
  * grid
  * gridExtra
  * ROCR
  * mgcv
  * ggmcmc


## Stan version for the sampling in the paper

The MCMC samples reported in the paper were all performed using Stan version 2.2.0.

## License

Copyright 2013, Thomas Wallis.

This is free, open source software released under the GPL-3 License. See LICENSE.txt for details.

Use this software at your own risk.
