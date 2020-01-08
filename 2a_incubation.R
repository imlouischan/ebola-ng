## incubation period ############################################################

# results from WHO, \cite{team2014ebola}
tauES_mean = 9.1; tauES_var = 7.3^2;
# shape & scale parameters, converted using continuous Gamma
thetaES = c(tauES_mean^2 / tauES_var, tauES_var / tauES_mean)

# discretization
pmf_ES = function(tauES) pgamma(tauES+1, shape = thetaES[1], scale = thetaES[2]) - pgamma(tauES, shape = thetaES[1], scale = thetaES[2])
cdf_ES = function(tauES) pgamma(tauES+1, shape = thetaES[1], scale = thetaES[2])

## plot ########################################################################

# plot name
file_plot <- paste0("figure/2a_incubation")
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 11, height = 8.5)

# observed incubation period
plot(prop.table(table(tauiES)), 
     xlim = c(0, max(tauiES, 14)), 
     cex.axis = 1.5, 
     xaxt = "n", 
     xlab = "", 
     ylab = "", 
     col = "black", lwd = 3)
axis(1, at = 0:max(tauiES, 14), labels = 0:max(tauiES, 14), cex.axis = 1.5)
title(xlab = "Incubation period", line = 2.5, cex.lab = 2)
title(ylab = "Probability", line = 2.5, cex.lab = 2)

# estimated pmf
points(0:max(tauiES, 14), pmf_ES(tauES = 0:max(tauiES, 14)), 
       cex = 2, 
       lwd = 3, 
       col = "blue")

# legend
legend("topright", 
       legend = c("data", "gamma"), 
       col = c("black", "blue"), 
       pch = c(NA, 1), 
       lty = c(1, 0), 
       lwd = 3, 
       cex = 2
)

# save as eps
dev.off()