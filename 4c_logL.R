## full log-likelihood (or posterior) ##########################################

function (thetaSS, R0, size, delta, epsilon, Data) {
  
  # call the function: likelihood of offspring number
  L_ri = dget("4a_L_ri.R")
  # call the function: likelihood of serial interval
  L_si = dget("4b_L_si.R")
  
  # preallocate
  L_ri_prop = rep(NA ,nrow(Data))
  L_si_prop = rep(NA ,nrow(Data))
  
  # each case
  for (i in 1:length(Data$i)) {
    
    # likelihood of offspring using proposed parameters
    L_ri_prop[i] = L_ri(i, thetaSS, R0, size, delta, epsilon, Data)
    # likelihood of serial intervals using proposed parameters
    L_si_prop[i] = L_si(i, thetaSS, R0, epsilon, Data)
    
  }
  
  # logL_r of all cases
  ( logL_r_prop = sum(log(L_ri_prop)) )
  
  # logL_s of all cases
  ( logL_s_prop = sum(log(L_si_prop)) )
  
  # log-likelihood
  ( logL_prop = logL_r_prop + logL_s_prop )
  
  # output
  if ( all(is.finite(logL_prop)) ) {
    return(logL_prop)
  } else { # uniform probability
    return(-1e10)
  }
  
}