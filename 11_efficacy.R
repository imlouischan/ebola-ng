## plot: efficacy ##############################################################

source("3b_pmf_s.R") # serial interval distribution

# hospitalization probability
ph = 1

# preallocate
type_data = rep(NA, nrow(Data))
type_col = rep(NA, nrow(Data))

# preallocate
efficacy = matrix(rep(NA, nrow(samples)*nrow(Data)), ncol = nrow(Data))

## calculation #################################################################
for (n in 1:nrow(samples)) { # 4000 samples
  
  # sampled parameters
  epsilon = samples$epsilon[n]
  thetaSS = c(samples$thetaSS1[n], samples$thetaSS2[n])
  
  for (i in 1:length(Data$i)) { # 20 cases
    
    # (1)(2)(3)(4) types
    if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
      type_data[i] = "(1) -OO-"
      type_col[i] = "chartreuse3"
      # data of observed cases
      tauSH = Data$tauiSH[i]
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = epsilon_h_i
      
    } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) ====
      type_data[i] = "(2) -OM-"
      type_col[i] = "cornflowerblue"
      # intervals
      a = 0
      b = Data_BH2(i) - Data$tiS[i]
      # change of variable
      tauSH = a:b
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = sum(epsilon_h_i*pmf_SH(tauSH))/sum(pmf_SH(tauSH))
      
    } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
      type_data[i] = "(3) OMO-"
      type_col[i] = "darkgoldenrod1"
      # intervals
      a = 0
      b = Data_BS2(i) - Data$tiE[i]
      # change of variable
      tauES = a:b
      tauSH = Data$tiH[i] - Data$tiE[i] - tauES
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = sum(epsilon_h_i*pmf_ES(tauES))/sum(pmf_ES(tauES))
      
    } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
      type_data[i] = "(4) MMO-"
      type_col[i] = "peachpuff3"
      # intervals
      a = Data$tiH[i] - Data_BS2(i)
      b = Data$tiH[i] - Data_BS1(i)
      # change of variable
      tauSH = a:b
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = sum(epsilon_h_i*pmf_SH(tauSH))/sum(pmf_SH(tauSH))
      
    }
  }
}

## plot ########################################################################
# plot name
file_plot <- paste0("figure/11_efficacy", 
                    "_N", N, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 11, height = 8.5)

# efficacy
boxplot(efficacy, 
        col = type_col, 
        lwd = 2, 
        cex.axis = 1.5, 
        xlab = "", 
        ylab = "",
        ylim = c(0, 1), yaxt = "n", 
        horizontal = TRUE)
# title(ylab = "Case ID (Infector ID)", line = 2.5, cex.lab = 2)
title(xlab = "Protective effect", line = 2.5, cex.lab = 2)

# Case ID & Infector ID
axis(2, at = c(1:nrow(Data)), las = 2, cex.axis = 1.0, 
     labels = paste(1:20, "(", Data$vi, ")", sep=""))

# effectiveness
abline(v = quantile(samples$epsilon, c(0.25, 0.5, 0.75)),
       lwd = c(2, 2, 2), lty = c("dashed", "solid", "dashed"), col = "red")
abline(v = c(min(samples$epsilon), max(samples$epsilon)),
       lwd = 2, lty = "dotted", col = "darkred")

# legend
legend("topright", 
       legend = c("(1)-OO-", "(2)-OM-", "(3)OMO-", "(4)MMO-"), 
       fill = c("chartreuse3", "cornflowerblue", "darkgoldenrod1", "peachpuff3"), 
       horiz = F, cex = 0.8, text.width = 0.063)

# save as eps
dev.off()