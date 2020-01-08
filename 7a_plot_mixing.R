## plot: mixing ################################################################
# plot name
file_plot <- paste0("figure/7a_mixing", 
                    "_N", N, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 11, height = 8.5)

# alternative models
samples_index = 4:11
# if (distr_r != "nbinom") { samples_index = setdiff(samples_index, 7) }
# if (!delta_TRUE)         { samples_index = setdiff(samples_index, 8) }

if (length(samples_index) == 6) {
  par(mfrow=c(3,2)) # subplot of 3x2
} else {
  par(mfrow=c(4,2)) # subplot of 4x2 
}

for (i in colnames(samples)[samples_index]) {
  plot(samples[, i], xlab = "", ylab = i) # plot mixing
}

par(mfrow=c(1,1)) # back to one plot

# save as eps
dev.off()