## model comparison ############################################################

i = 0
comp_model = data.frame(
  delta.TRUE   = rep(NA, 1),
  distr.r      = rep(NA, 1),
  distr.s      = rep(NA, 1),
  k        = rep(NA, 1),
  logMLE   = rep(NA, 1),
  AIC      = rep(NA, 1),
  AIC_diff = rep(NA, 1),
  BIC      = rep(NA, 1),
  BIC_diff = rep(NA, 1),
  ML       = rep(NA, 1),
  stringsAsFactors = FALSE # undefine levels
)

delta_TRUE <- 1
for (distr_r in c("geom", "pois", "nbinom")) { # alternative offspring distributions
  for (distr_s in c("dgamma", "dweibull")) { # alternative serial interval distributions
    
    # samples of each alternative model
    samples = samples_alt[samples_alt$delta.TRUE   == delta_TRUE &
                            samples_alt$distr.r    == distr_r &
                            samples_alt$distr.s    == distr_s, ]
    
    ## AIC #########################################################################
    
    k = 5 # number of parameters
    if (distr_r == "nbinom") { k = k + 1 }
    if (delta_TRUE == 1    ) { k = k + 1 }
    
    # MLE
    logMLE = max(samples$logL)
    
    AIC = 2*k - 2*logMLE
    
    ## BIC #########################################################################
    
    n = nrow(Data) # sample size
    
    BIC = log(n)*k - 2*logMLE
    
    ## PHM estimator ###############################################################
    
    # likelihood
    L = exp(samples$logL)
    
    # marginal likelihood
    ML = 1/sum(1/L) / length(L)
    
    # output
    i = i + 1
    comp_model[i, ] = list(delta_TRUE, distr_r, distr_s, k, logMLE, AIC, NA, BIC, NA, ML)
    
  }
} 

# differences in AIC & BIC
comp_model$AIC_diff = comp_model$AIC - min(comp_model$AIC)
comp_model$BIC_diff = comp_model$BIC - min(comp_model$BIC)

# model posterior
comp_model$pM = comp_model$ML / sum(comp_model$ML) * 100

# output
print( format(comp_model, digits = 1, nsmall = 1) )