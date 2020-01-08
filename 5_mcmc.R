## MCMC ########################################################################

function (N, chain, delta_TRUE) {
  
  # control randomness, choose whatever you like
  set.seed(chain + 2018)
  
  # call the function: log-posterior (or log-likelihood)
  logL = dget("4c_logL.R")
  # call the function: transmission probability
  pij = dget("4d_pij.R")
  
  # preallocate
  s.thetaSS = matrix(rep(c(NA, NA) ,N), ncol = 2)
  s.R0      = rep(NA ,N)
  s.size    = rep(NA ,N)
  s.delta   = rep(NA ,N)
  s.epsilon = rep(NA ,N)
  s.vi      = rep(NA ,N)
  s.logL    = rep(NA ,N)
  
  # counter of acceptance
  count_acc = 0
  # timer
  timer_start = proc.time()
  
  # subchain number
  N.AM = N/100
  
  ## start MCMC ##################################################################
  for (n in 1:N) {
    
    if ( (n-1)%%N.AM == 0 & n > N.AM & n-1 <= N/2 ) { # update proposal distribution, during burn-in
      # index of previous samples
      n.AM = (n - N.AM):(n - 1)
      # update proposal using previous samples
      thetaSS_sd = rep(NA ,2)
      thetaSS_sd[1] = sd(s.thetaSS[n.AM, 1]) + 1e-6
      thetaSS_sd[2] = sd(s.thetaSS[n.AM, 2]) + 1e-6
      R0_sd         = sd(s.R0[n.AM])         + 1e-6
      size_sd       = sd(s.size[n.AM])       + 1e-6
      delta_sd      = sd(s.delta[n.AM])      + 1e-6
      epsilon_sd    = sd(s.epsilon[n.AM])    + 1e-6
    }
    
    if (n == 1) { # initial parameters
      
      # proposing parameters (NOT from uniform prior)
      thetaSS = rep(NA ,2)
      thetaSS[1] = runif(1, 0, 20)
      thetaSS[2] = runif(1, 0, 10)
      R0         = runif(1, 0, 19)
      size       = runif(1, 0, 5)
      delta      = runif(1, 0, 0.5)
      epsilon    = runif(1, 0, 1)
      
    } else if ( n <= N.AM ) { # proposal distribution, at the very begining
      
      # Gaussian proposal distribution
      thetaSS = rnorm(2, mean = s.thetaSS[n - 1, ], sd = c(5, 1))
      R0      = rnorm(1, mean = s.R0[n - 1],        sd = 1)
      size    = rnorm(1, mean = s.size[n - 1],      sd = 0.5)
      delta   = rnorm(1, mean = s.delta[n - 1],     sd = 0.02)
      epsilon = rnorm(1, mean = s.epsilon[n - 1],   sd = 0.3)
      
    } else { # proposal distribution
      
      # proposal distribution using previous samples
      thetaSS = rnorm(2, mean = s.thetaSS[n - 1, ], sd = thetaSS_sd)
      R0      = rnorm(1, mean = s.R0[n - 1],        sd = R0_sd)
      size    = rnorm(1, mean = s.size[n - 1],      sd = size_sd)
      delta   = rnorm(1, mean = s.delta[n - 1],     sd = delta_sd)
      epsilon = rnorm(1, mean = s.epsilon[n - 1],   sd = epsilon_sd)
      
    }
    
    # uniform prior, for ddgamma & ddweibull & dnbinom
    if (thetaSS[1] > 0 & thetaSS[1] < 1000 & 
        thetaSS[2] > 0 & thetaSS[2] < 100 & 
        R0         > 0 & R0         < 19  & 
        size       > 0 & size       < 100 & 
        delta      > 0 & delta      < 1   & 
        epsilon    > 0 & epsilon    < 1   ){
      
      # alternative models
      if (!delta_TRUE) delta = 0
      if (distr_r != "nbinom") size = 1
      
      ## proposal distribution of missing infector with transmission prob $p_{ij}$ ####
      
      # contacts with Case 20
      wi = c(5, 6, 7, 10, 11, 12)
      # proposed infector
      ( Data$vi[20] = as.integer(sample(wi, size = 1, prob = pij(thetaSS, R0, epsilon, Data))) )
      
      ## likelihood calculation ######################################################
      
      # log-likelihood
      logL_prop = logL(thetaSS, R0, size, delta, epsilon, Data)
      
      # make sure initial parameters do not produce "logL = -Inf"
      if (n == 1) {
        while ( !is.finite(logL_prop) ) {
          # proposing parameters (not from uniform prior)
          thetaSS = rep(NA ,2)
          thetaSS[1] = runif(1, 0, 20)
          thetaSS[2] = runif(1, 0, 10)
          R0         = runif(1, 0, 19)
          size       = runif(1, 0, 5)
          delta      = runif(1, 0, 0.5)
          epsilon    = runif(1, 0, 1)
          # log-likelihood
          logL_prop = logL(thetaSS, R0, size, delta, epsilon, Data)
        }  
      } 
      
    } else { # prior, out of range
      logL_prop = -Inf
    }
    
    ## accept or reject ############################################################
    
    # acceptance probability
    ( logA = min(0, logL_prop - s.logL[n - 1]) )
    ( logu = log(runif(1)) )
    
    # accept or reject
    if (n == 1 | logu < logA) { # initial value and accept
      s.thetaSS[n, ] = thetaSS
      s.R0[n]        = R0
      s.size[n]      = size
      s.delta[n]     = delta
      s.epsilon[n]   = epsilon
      s.vi[n]        = Data$vi[20]
      s.logL[n]      = logL_prop
      # counter of acceptance
      count_acc = count_acc + 1
    } else { # reject
      s.thetaSS[n, ] = s.thetaSS[n - 1, ]
      s.R0[n]        = s.R0[n - 1]
      s.size[n]      = s.size[n - 1]
      s.delta[n]     = s.delta[n - 1]
      s.epsilon[n]   = s.epsilon[n - 1]
      s.vi[n]        = s.vi[n - 1]
      s.logL[n]      = s.logL[n - 1]
    }
    
    # counter
    if ( (n*10)%%N == 0 ) {
      percent = paste(n, "/", N, " = ", n*100/N, "%", sep="")
      timer = paste(round((proc.time() - timer_start)/60), "min", sep="")
      acc_rate = paste("acc.rate = ", round(count_acc*1000/N), "%", sep=""); count_acc = 0
      print(c(percent, timer[1], acc_rate))
    }
  }
  
  ## thinning and save ###########################################################
  
  # output sample number
  N.thin = 1000
  # thinning & burn-in
  thin = seq(N/2 + 1, N, by = N/N.thin/2)
  
  # data frame
  samples = data.frame(thetaSS1 = s.thetaSS[thin, 1], 
                       thetaSS2 = s.thetaSS[thin, 2], 
                       R0       = s.R0[thin], 
                       size     = s.size[thin], 
                       delta    = s.delta[thin], 
                       epsilon  = s.epsilon[thin], 
                       vi       = s.vi[thin], 
                       logL     = s.logL[thin])
  
  # save samples
  file_samples = paste("samples", 
                       "_N", N, 
                       "_d", delta_TRUE, 
                       "_R", distr_r, 
                       "_S", distr_s,
                       "_C", chain, 
                       ".Rdata", sep="")
  save(samples, file = file_samples)
  
  # output
  return(samples)
  
}