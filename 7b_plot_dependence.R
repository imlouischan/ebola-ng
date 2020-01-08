## plot: dependence ############################################################
# plot name
file_plot <- paste0("figure/7b_pairwise", 
                    "_N", N, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# alternative models
samples_index = 4:9
if (distr_r != "nbinom") { samples_index = setdiff(samples_index, 7) }
if (!delta_TRUE)         { samples_index = setdiff(samples_index, 8) }

# plot parameter dependence
library(GGally)
print( 
  GGally::ggpairs(samples[, samples_index], 
                  lower = list(continuous = GGally::wrap("points", alpha = 0.1)), 
                  upper = list(continuous = "cor"),
                  diag = list(continuous = "densityDiag")) + 
    ggplot2::theme_bw()
)

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, width = 11, height = 8.5, device = cairo_ps)