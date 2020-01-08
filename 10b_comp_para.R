## parameter comparison ########################################################
# plot name
file_plot <- paste0("figure/10b_sensitivity", 
                    "_N", N, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# samples of alternative models
for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
  samples_delta = samples_alt[samples_alt$distr.r    == distr_r &
                                samples_alt$distr.s    == distr_s, ]
}

# plots
library(ggplot2)

p1 = ggplot(samples_delta, aes(x = thetaSS1, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p2 = ggplot(samples_delta, aes(x = thetaSS2, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p3 = ggplot(samples_delta, aes(x = R0, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p4 = ggplot(samples_delta, aes(x = size, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p5 = ggplot(samples_delta[samples_delta$delta.TRUE == 1, ], aes(x = delta, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p6 = ggplot(samples_delta, aes(x = epsilon, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p7 = ggplot(samples_delta, aes(x = reorder(vi, vi, function(x) - length(x)), linetype = delta.TRUE)) +
  geom_bar(aes(y = (..count..)/(sum(..count..)/2)), position = "dodge", col = "black") +
  xlab("Infector ID") + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  scale_linetype_discrete(name = "", labels = c("Variable Rt", "Constant R0")) + 
  theme(legend.position = c(0.75, 0.7), legend.title = element_blank())

p8 = ggplot(samples_delta, aes(x = logL, linetype = delta.TRUE)) + 
  geom_density() + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

# multiple plots
print( # ggplot requires an explicit print
  cowplot::plot_grid( p1, p2, p3,     p5, p6, p7,       nrow = 3, 
                      # p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4, 
                      # labels = "AUTO", 
                      align = "hv"
  )
)

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, width = 11, height = 8.5)