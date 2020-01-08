## likelihood of offspring number ##############################################

function (i, thetaSS, R0, size, delta, epsilon, Data) {
  
  # preallocate
  L_ri_1 = 0
  L_ri_2 = 0
  
  # observed offspring number
  ri_df = as.data.frame(table(factor(Data$vi, levels = Data$i)))
  ( ri = ri_df$Freq )
  
  # hospitalization probability
  ph = 1
  
  # categories in L_r (max: 7 types)
  # i: "xOOx", "xOMx", "OMOx", "MMOx", "OMMO", "OMMM", "MMMO", ("MMMM")
  if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
    type.i = "(1) xOOx"
    # data of observed cases
    tauSH = Data$tauiSH[i]
    # exponentially decreasing R
    R = R0*exp(-delta*Data$tiS[i])
    # L_r of each case
    L_ri_1 = pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
    L_ri_2 = 1
  } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) ====
    type.i = "(2) xOMx"
    # intervals
    a = 0
    b = Data_BH2(i) - Data$tiS[i]
    # change of variable
    tauSH = a:b
    # exponentially decreasing R
    R = R0*exp(-delta*Data$tiS[i])
    # numerator and denominator
    L_ri_1 = sum( pmf_SH(tauSH)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size) )
    L_ri_2 = sum( pmf_SH(tauSH) )
  } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
    type.i = "(3) OMOx"
    # intervals
    a = 0
    b = Data_BS2(i) - Data$tiE[i]
    # change of variable
    tauES = a:b
    tauSH = Data$tiH[i] - Data$tiE[i] - tauES
    # exponentially decreasing R
    R = R0*exp(-delta*(Data$tiE[i] + tauES))
    # numerator and denominator
    L_ri_1 = sum( pmf_ES(tauES)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size) )
    L_ri_2 = sum( pmf_ES(tauES) )
  } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
    type.i = "(4) MMOx"
    # intervals
    a = Data$tiH[i] - Data_BS2(i)
    b = Data$tiH[i] - Data_BS1(i)
    # change of variable
    tauSH = a:b
    # exponentially decreasing R
    R = R0*exp(-delta*(Data$tiH[i] - tauSH))
    # numerator and denominator
    L_ri_1 = sum( pmf_SH(tauSH)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size) )
    L_ri_2 = sum( pmf_SH(tauSH) )
  } else { print(i); break } # error(?)
  
  # the follow categories do not exist in Nigeria data, \cite{folarin2016ebola}
  # if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) {print(c(i,"OMMO",5))};
  # if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) {print(c(i,"OMMM",6))};
  # if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) {print(c(i,"MMMO",7))};
  # if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) {print(c(i,"MMMM",8))};
  
  # resulting likelihood
  ( L_ri = L_ri_1/L_ri_2 )
  
  # output
  return(L_ri)
  
}