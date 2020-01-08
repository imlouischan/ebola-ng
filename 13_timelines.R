## plot: timelines #############################################################
# http://www.statmethods.net/graphs/line.html
# http://www.statmethods.net/advgraphs/parameters.html
# http://www.statmethods.net/graphs/dot.html

# plot name
file_plot <- paste0("figure/13_timelines")
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 11, height = 8.5)

# empty plot
plot(0, 0, type = "n",
     xlim = c(0, Data_tEND), 
     ylim = c(1, nrow(Data)), 
     xlab = "", 
     ylab = "", 
     xaxt = "n", 
     yaxt = "n") 
# title(ylab = "Case ID (Infector ID)", line = 2.5, cex.lab = 2)
title(xlab = "Day", line = 2.5, cex.lab = 2)
axis(1, at = seq(0, Data_tEND, 10), 
     labels = seq(0, Data_tEND, 10), 
     cex.axis = 1.5)

# Case ID & Infector ID
axis(2, at = c(1:nrow(Data)), las = 2, cex.axis = 1.0, 
     labels = paste(Data$i, "(", Data$vi, ")", sep=""))

# the enlarged size of estimated time points (x: cross) for better visualization
enlarge = 10

for (i in 1:nrow(Data)) { # 20 cases
  
  # background horizontal lines
  lines(c(c(0, Data_tEND)), c(i, i),
        lwd = 1, lty = "dotted", col = "gray")
  
  # all observed time points
  points(Data$tiE[i], i,
         col = "darkgray", lwd = 2, pch = 1)
  points(Data$tiS[i], i,
         col = "red"     , lwd = 2, pch = 1)
  points(Data$tiH[i], i,
         col = "blue"    , lwd = 2, pch = 1)
  points(Data$tiD[i], i,
         col = "black"   , lwd = 2, pch = 1)
  
  # estimated time points, with PMFs
  if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
    type.i = "(1) xOOx"
    
  } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) ====
    type.i = "(2) xOMx"
    # intervals
    a = 0
    b = Data_BH2(i) - Data$tiS[i]
    # change of variable
    tauSH = a:b
    
    # possible time points
    points(Data$tiS[i] + tauSH, rep(i, length(tauSH)),
           col = "blue", lwd = 2, pch = 4,
           cex = pmf_SH(tauSH)/sum(pmf_SH(tauSH)) * enlarge + 1e-20) # avoid zero
    
  } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
    type.i = "(3) OMOx"
    # intervals
    a = 0
    b = Data_BS2(i) - Data$tiE[i]
    # change of variable
    tauES = a:b
    tauSH = Data$tiH[i] - Data$tiE[i] - tauES
    
    # possible time points
    points(Data$tiE[i] + tauES, rep(i, length(tauES)),
           col = "red", lwd = 2, pch = 4,
           cex = pmf_ES(tauES)/sum(pmf_ES(tauES)) * enlarge)
    
  } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
    type.i = "(4) MMOx"
    # intervals
    a = Data$tiH[i] - Data_BS2(i)
    b = Data$tiH[i] - Data_BS1(i)
    # change of variable
    tauSH = a:b
    
    # possible time points
    points(Data$tiH[i] - tauSH, rep(i, length(tauSH)),
           col = "red", lwd = 2, pch = 4,
           cex = pmf_SH(tauSH)/sum(pmf_SH(tauSH)) * enlarge)
  }
}

# legend, both observed and estimated
legend("bottomright", legend = c("Date of exposure", 
                                 "Date of illness onset", 
                                 "Date of illness onset (estimated)", 
                                 "Date of hospitalization", 
                                 "Date of hospitalization (estimated)", 
                                 "Date of death"),
       col = c("darkgray", "red", "red", "blue", "blue", "black"), 
       pch = c(1, 1, 4, 1, 4, 1), 
       lwd = 2, 
       lty = 0, 
       cex = 1.5
)

# save as eps
dev.off()