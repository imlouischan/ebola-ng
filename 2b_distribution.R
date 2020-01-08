## MLE parameter estimation of time lag distributions ##########################

## time lag
## (1) tauSH: time lags from illness onset to hospitalization
## (2) tauSS: serial intervals (in the presence of isolation)

## alternative distributions (pmf) and MLE
## (a) gamma distribution
## (b) weibull distribution

# alternative distributions
pmf_ddgamma   = function(tau, theta)   pgamma(tau+1, shape = theta[1], scale = theta[2]) -   pgamma(tau, shape = theta[1], scale = theta[2])
cdf_ddgamma   = function(tau, theta)   pgamma(tau+1, shape = theta[1], scale = theta[2])
pmf_ddweibull = function(tau, theta) pweibull(tau+1, shape = theta[1], scale = theta[2]) - pweibull(tau, shape = theta[1], scale = theta[2])
cdf_ddweibull = function(tau, theta) pweibull(tau+1, shape = theta[1], scale = theta[2])

pmf_alt = function(distr_index) {
  if (distr_index == 1) return(pmf_ddgamma)
  if (distr_index == 2) return(pmf_ddweibull)
}

# alternative observed time lags
taui_alt = function(taui_index) {
  if (taui_index == 1) return(tauiSH)
  if (taui_index == 2) return(tauiSS)
}
taui_name <- c("Time from illness onset to hospitalization", "Serial interval")

# preallocate
MLE_alt   = list()
distr_MLE = list()
theta_MLE = list()

for (taui_index in 1:2) { # different observed time lags
  
  # plot name
  file_plot <- paste0("figure/2b_distribution", 
                      "_", taui_index)
  # save as eps
  setEPS()
  print( name_file <- paste0(file_plot, ".eps") )
  postscript(name_file, width = 11, height = 8.5)
  
  # preallocate
  MLE_alt[[taui_index]] = data.frame(
    dgamma   = rep(NA, 3),
    dweibull = rep(NA, 3),
    row.names = c("nlogL", "theta1", "theta2")
  )
  
  # plot observed time lags
  plot(prop.table(table(taui_alt(taui_index))),
       xlim = c(0, max(taui_alt(taui_index), 14)),
       cex.axis = 1.5, 
       xaxt = "n", 
       xlab = "", 
       ylab = "", 
       col = "black", lwd = 3)
  axis(1, at = 0:max(taui_alt(taui_index), 14), labels = 0:max(taui_alt(taui_index), 14), cex.axis = 1.5)
  title(xlab = taui_name[taui_index], line = 2.5, cex.lab = 2)
  title(ylab = "Probability", line = 2.5, cex.lab = 2)
  
  for (distr_index in 1:2) { # alternative distributions
    
    # negative log-likelihood
    nlogL = function(theta) -sum(log(pmf_alt(distr_index)(taui_alt(taui_index), theta)))
    
    # parameter estimation using optim
    theta_optim = suppressWarnings( optim(c(1, 1), nlogL) )
    MLE_alt[[taui_index]]["nlogL", distr_index] = theta_optim$value
    MLE_alt[[taui_index]][c("theta1", "theta2"), distr_index] = theta_optim$par
    
    # plot estimated pmf
    points(0:max(taui_alt(taui_index), 14), 
           pmf_alt(distr_index)(tau = 0:max(taui_alt(taui_index), 14), theta = theta_optim$par), 
           cex = 2, 
           lwd = 3, 
           pch = distr_index, 
           col = "gray")
  }  
  
  # checking
  print(MLE_alt[[taui_index]])
  
  # the best distribution and estimated parameters
  distr_MLE[[taui_index]] = which.min(MLE_alt[[taui_index]]["nlogL", ])
  theta_MLE[[taui_index]] = MLE_alt[[taui_index]][c("theta1", "theta2"), distr_MLE[[taui_index]]]
  
  # best pmf
  if (taui_index == 1) pmf_SH = function(tau) pmf_alt(distr_MLE[[1]])(tau, theta_MLE[[1]]) 
  if (taui_index == 2) pmf_SS = function(tau) pmf_alt(distr_MLE[[2]])(tau, theta_MLE[[2]])
  
  # plot best pmf
  points(0:max(taui_alt(taui_index), 14), 
         pmf_alt(distr_MLE[[taui_index]])(tau = 0:max(taui_alt(taui_index), 14), theta = theta_MLE[[taui_index]]), 
         cex = 2, 
         lwd = 3, 
         pch = distr_MLE[[taui_index]], 
         col = "red")
  
  # legend
  ifelse(distr_MLE[[taui_index]]==1, 
         # gamma is the best
         legend("topleft", 
                legend = c("data", "gamma (best fit)", "Weibull"), 
                col = c("black", "red", "gray"), 
                pch = c(NA, 1, 2), 
                lty = c(1, 0, 0), 
                lwd = 3, 
                cex = 2
         ), 
         # Weibull is the best
         legend("topright", 
                legend = c("data", "gamma", "Weibull (best fit)"), 
                col = c("black", "gray", "red"), 
                pch = c(NA, 1, 2), 
                lty = c(1, 0, 0), 
                lwd = 3, 
                cex = 2
         ))
  
  # save as eps
  dev.off()
  
} # different observed time lags