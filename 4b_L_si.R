## likelihood of serial interval ###############################################

function (i, thetaSS, R0, epsilon, Data) {
  
  # preallocate
  L_si_1 = 0
  L_si_2 = 0
  
  # infector of i
  vi = Data$vi[i]
  
  # hospitalization probability
  ph = 1
  
  # categories in L_s (max: 28 types)
  # vi: "xOOx", "xOMx", "OMOx", "MMOx", "OMMO", "OMMM", "MMMO"
  #  i: "xOxx", "OMxx", "MMOx", "MMMO"
  if (vi == 0) { # the index case / imported cases
    type.vi = "(0)"
    type.i = "(0)"
    # L_s of each case
    L_si_1 = 1
    L_si_2 = 1
  } else { # except the index case / imported cases
    if (                      Data_obs[vi,2] &  Data_obs[vi,3]                     ) { # vi: xOOx (1) ====
      type.vi = "(1) xOOx"
      if (             Data_obs[i,2]                        ) { # i: xOxx (I) ----
        type.i = "(I) xOxx"
        # data of observed cases
        tauSS = Data$tiS[i] - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # L_s of each case
        L_si_1 = pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
        L_si_2 = 1
      }
      else if ( Data_obs[i,1] & !Data_obs[i,2]                        ) { # i: OMxx (II) ----
        type.i = "(II) OMxx"
        # intervals
        a = 0
        b = Data_BS2(i) - Data$tiE[i]
        # change of variable
        tauES = a:b
        tauSS = Data$tiE[i] + tauES - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # numerator and denominator
        L_si_1 = sum( pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_ES(tauES) )
      }
      else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (III) ----
        type.i = "(III) MMOx"
        # intervals
        a = Data$tiH[i] - Data_BS2(i)
        b = Data$tiH[i] - Data$tiS[vi]
        # change of variable
        tauiSH = a:b
        tauSS = Data$tiH[i] - tauiSH - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # numerator and denominator
        L_si_1 = sum( pmf_SH(tauiSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_SH(tauiSH) )
      }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO (IV)
        # not exist in Nigeria data 
      }
    }
    else if (                      Data_obs[vi,2] & !Data_obs[vi,3]                     ) { # vi: xOMx (2) ====
      type.vi = "(2) xOMx"
      if (             Data_obs[i,2]                        ) { # i: xOxx (I) ----
        type.i = "(I) xOxx"
        # intervals
        a = 0
        b = Data_BH2(vi) - Data$tiS[vi]
        # change of variable
        tauSS = Data$tiS[i] - Data$tiS[vi]
        tauSH = a:b
        # numerator and denominator
        L_si_1 = sum( pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_SH(tauSH) )
      }
    }
    else if ( Data_obs[vi,1] & !Data_obs[vi,2] &  Data_obs[vi,3]                     ) { # vi: OMOx (3) ====
      type.vi = "(3) OMOx"
      if (             Data_obs[i,2]                        ) { # i: xOxx (I) ----
        type.i = "(I) xOxx"
        # intervals
        a = 0
        b = Data_BS2(vi) - Data$tiE[vi]
        # change of variable
        tauES = a:b
        tauSS = Data$tiS[i] - Data$tiE[vi] - tauES
        tauSH = Data$tiH[vi] - Data$tiE[vi] - tauES
        # numerator and denominator
        L_si_1 = sum( pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_ES(tauES) )
      }
    }
    else if (!Data_obs[vi,1] & !Data_obs[vi,2] &  Data_obs[vi,3]                     ) { # vi: MMOx (4) ====
      type.vi = "(4) MMOx"
      if (             Data_obs[i,2]                        ) { # i: xOxx (I) ----
        type.i = "(I) xOxx"
        # intervals
        a = Data$tiH[vi] - Data_BS2(vi)
        b = Data$tiH[vi] - Data_BS1(vi)
        # change of variable
        tauSS = Data$tiS[i] - Data$tiH[vi] + (a:b)
        tauSH = a:b
        # numerator and denominator
        L_si_1 = sum( pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_SH(tauSH) )
      }
    }
    # the follow categories do not exist in Nigeria data, \cite{folarin2016ebola}
    # vi: "OMMO", "OMMM", "MMMO"
    # if ( Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] &  Data_obs[vi,4]) { 
    # if ( Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] & !Data_obs[vi,4]) { 
    # if (!Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] &  Data_obs[vi,4]) {     
    # if (!Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] & !Data_obs[vi,4]) { 
  }
  
  # resulting likelihood
  ( L_si = L_si_1/L_si_2 )
  
  # output each case
  return(L_si)
  
}