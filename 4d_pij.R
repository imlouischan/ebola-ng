## transmission probability $p_{ij}$ ###########################################

function (thetaSS, R0, epsilon, Data) {
  
  # call the function: likelihood of serial interval
  L_si = dget("4b_L_si.R")
  
  # contacts with Case 20
  w20 = c(5, 6, 7, 10, 11, 12)
  
  # preallocate
  s_hat20 = rep(NA ,length(w20))
  
  for (j in 1:length(w20)){
    
    # the missing infector
    Data$vi[20] = w20[j]
    
    # likelihood of serial interval using proposed parameters
    # actually the serial interval pmf
    s_hat20[j] = L_si(i = 20, thetaSS, R0, epsilon, Data)
    
  }
  
  # transmission probability
  p_20_j = s_hat20/sum(s_hat20)
  
  # output the transmission probability
  if ( all(is.finite(p_20_j)) ) {
    return(p_20_j)
  } else { # uniform probability
    return(NULL)
  }
  
}