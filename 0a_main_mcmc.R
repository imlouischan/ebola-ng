## new script ##################################################################
cat("\014") # clear the console
rm(list=ls()) # remove all variables
graphics.off() # close all plots
setwd("../ebola_ng/")

## main: MCMC ##################################################################

source("1a_data.R") # data set
source("1b_bound.R") # boundaries of missing time points

source("2a_incubation.R") # incubation period
source("2b_distribution.R") # estimation of time lag distributions

# alternative models
for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
  for (distr_r in c("geom", "pois", "nbinom")) { # alternative offspring distributions
    for (distr_s in c("dgamma", "dweibull")) { # alternative serial interval distributions
      
      # sample number of a chain
      N = 1e6 # can try the minimun N = 2e3
      
      # counter
      print(paste0("******************************", 
                   "_N", N, 
                   "_d", delta_TRUE, 
                   "_R", distr_r, 
                   "_S", distr_s, 
                   "******************************"))
      
      source("3a_pmf_r.R") # offspring distribution
      source("3b_pmf_s.R") # serial interval distribution
      
      # call the function to run one mcmc chain
      run_chain = dget("5_mcmc.R")
      run_chains = function(chain) {run_chain(N, chain, delta_TRUE)}
      
      # option A: local running
      # samples = run_chains(0)
      
      # option B: parallel running
      source("6a_parallel.R")
      
    }
  }
}

## main: summary ###############################################################

# preallocate
samples_alt = data.frame(NULL)

for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
  for (distr_r in c("geom", "pois", "nbinom")) { # alternative offspring distributions
    for (distr_s in c("dgamma", "dweibull")) { # alternative serial interval distributions
      
      
      # combine chains of each alternative model
      comb_chain = dget("6b_combine_chains.R")
      samples = comb_chain(N, 1:no_cores)
      
      # combine alternative models
      samples_alt = rbind(samples_alt, samples)
      
    }
  }
}

# save all combined samples
file_samples = paste0("samples_alt", 
                     "_N", N, 
                     ".Rdata")
save(samples_alt, file = file_samples)