## new script ##################################################################
cat("\014") # clear the console
rm(list=ls()) # remove all variables
graphics.off() # close all plots
setwd("../ebola_ng/")

## call data and results #######################################################

source("1a_data.R") # data set
source("1b_bound.R") # boundaries of missing time points

source("2a_incubation.R") # incubation period
source("2b_distribution.R") # estimation of time lag distributions

source("13_timelines.R") # plot timelines

# sample number of a chain
N = 1e6 # can try the minimun N = 2e3

# load all samples
file_samples = paste0("samples_alt", 
                     "_N", N, 
                     ".Rdata")
load(file = file_samples)

## plot (best model) ###########################################################

# the best model
delta_TRUE = 1; distr_r = "geom"; distr_s = "dgamma"; 

# mcmc samples of an alternative model
samples = samples_alt[samples_alt$delta.TRUE   == delta_TRUE &
                        samples_alt$distr.r    == distr_r &
                        samples_alt$distr.s    == distr_s, ]

source("7a_plot_mixing.R") # plot: mixing
source("7b_plot_dependence.R", echo = F) # plot: dependence

source("8a_plot_pmf_r.R") # plot: offspring distributions
source("8b_plot_pmf_s.R") # plot: serial interval distributions

source("9a_cal_pmf_r.R") # calculate: mean & SD of offspring distribution
source("9b_cal_pmf_s.R") # calculate: mean & SD of serial interval

source("11_efficacy.R") # plot: efficacy
source("12_infector.R", echo = F) # plot: possible infectors

## plot (alternative models) ###################################################

source("10a_comp_model.R") # plot: model comparison
source("10b_comp_para.R") # plot: parameter comparison