## combine chains ##############################################################

function (N, chains) { 
  
  # preallocate
  samples.all = data.frame(NULL)
  
  for (chain in chains) {
    
    # load samples
    file_samples = paste0("samples", 
                         "_N", N, 
                         "_d", delta_TRUE, 
                         "_R", distr_r, 
                         "_S", distr_s,
                         "_C", chain, 
                         ".Rdata")
    load(file = file_samples)
    
    # alternative model
    samples = cbind(
      delta.TRUE   = delta_TRUE, 
      distr.r      = distr_r,
      distr.s      = distr_s,
      samples)
    
    # factor instead of numeric
    samples$vi           = as.factor(samples$vi)
    samples$delta.TRUE   = as.factor(samples$delta.TRUE)
    
    # combine chains
    samples.all = rbind(samples.all, samples)
    
  }
  # output samples of all chains
  return(samples.all) 
}