############################################################################
# Bioramble
# A closer look at the fisherman's dilemma
# by Jesse Lipp
# Aug 7, 2015
############################################################################

# --------------------------------------------------------------------------
# Set up environment
# --------------------------------------------------------------------------
# clean-up
rm(list = ls())

# load libraries (install if needed)
if (!require(ggplot2)) {
  install.packages("ggplot2")
}
if (!require(reshape2)) {
  install.packages("reshape2")
}

# --------------------------------------------------------------------------
# prepare data
# --------------------------------------------------------------------------
# function implementing formula to calculate positive predictive value (PPV)
calc_ppv <- function(prior, power, alpha) {
  power * prior / (power * prior + alpha * (1 - prior))
}
# establish domain to calculate PPV
prior <- seq(0, 1, length.out = 100)
power <- seq(0, 1, length.out = 100)
# four commonly used alphas in high-throughput screening to calculate PPV contour
alpha <- c(0.1, 0.05, 0.01, 0.001)
# set up grid to calculate all combinations of prior and power
grid <- expand.grid(list(prior = prior, power = power))
# calculate PPV for all combinations of prior and power for all levels of alpha
ppv <- mapply(calc_ppv, grid$prior,  grid$power, MoreArgs = list(alpha = alpha))
# correct 0-division
ppv[is.na(ppv)] <- 0
# combine data into single data frame
rownames(ppv) <- alpha
grid <- cbind(grid, t(ppv))
# create long form for plotting with ggplot2
grid <- melt(grid,
             id = c("prior", "power"),
             value.name = "PPV",
             variable.name = "alpha")

# plot PPV contours as facets
ggplot(grid, aes(x = prior, y = power, fill = PPV)) +
  # color gradient from green to red
  geom_tile() + 
  scale_fill_gradient(low = "red", high = "green") +
  # contour line at PPV = 0.5
  stat_contour(aes(z = PPV), breaks = 0.5) +
  # facets
  facet_wrap(~ alpha, ncol = 2) +
  # rectangle specifying the region of most high-throughput screens
  geom_rect(xmin = 0, xmax = 0.25, ymin = 0, ymax = 0.25, fill = NA, color = "black", linetype = 3) + 
  annotate("text", label = "HTS", x = 0.125, y = 0.125, size = 10) +
  # labeling
  xlab("Prior") + ylab("Power") +
  # theme tweaking
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black"))
