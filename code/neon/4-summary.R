# 4-summary.R: this script summarizes results from the Central Florida 
#              bird case study.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)
library(MCMCvis)
library(tidyverse)
library(sf)
library(stars)
library(viridis)
library(cowplot)
library(ggspatial)
library(unmarked)

# Reorder neonDWP data following order used to fit model ------------------
start.sp <- c('EATO', 'BACS')
sp.names <- dimnames(neonDWP$y)[[1]]
# Other species codes
indices <- rep(NA, length(start.sp))
for (i in 1:length(indices)) {
  indices[i] <- which(sp.names == start.sp[i])
}
indices.other <- 1:nrow(neonDWP$y)
indices.other <- indices.other[-indices]
# Create the ordered y data frame
neonDWP$y <- neonDWP$y[c(indices, indices.other), , ]
# Updated species codes
sp.names <- sp.names[c(indices, indices.other)]
dimnames(neonDWP$y)[[1]] <- sp.names

# Read in results for comparison ------------------------------------------
# The multi-species results files are too large for GitHub, so condensed
# versions of those objects are read in afterwards. These objects can be 
# obtained by running the "main" scripts, or contact Jeff Doser (doserjef@msu.edu)
# load("results/neon-msDS-results.rda")
# out.msDS <- out
# load("results/neon-lfMsDS-results.rda")
# out.lfMsDS <- out
# load("results/neon-sfMsDS-results.rda")
# out.sfMsDS <- out
# 
# # Quick summaries ---------------------------------------------------------
# summary(out.msDS, level = 'community')
# summary(out.lfMsDS, level = 'community')
# summary(out.sfMsDS, level = 'community')

# Read in the smaller model objects available on GitHub -------------------
# Loads the following objects for each model type (* is replaced by the function name): 
#   - alpha.comm.samples.*: community-level detection effects
#   - alpha.samples.*: species-level detection effects
#   - beta.comm.samples.*: community-level abundance effects
#   - beta.samples.*: species-level abundance effects
#   - lambda.samples.*: factor loadings
#   - rhat.*: potential scale reduction factors for model parameters
#   - tau.sq.alpha.samples.*: community-level detection variances
#   - tau.sq.beta.samples.*: community-level abundance variances
# Look at convergence of all model objects
load("results/neon-msDS-small-results.rda")
rhat.msDS
load("results/neon-lfMsDS-small-results.rda")
rhat.lfMsDS
load("results/neon-sfMsDS-small-results.rda")
rhat.sfMsDS

# Calculate WAIC (takes a few minutes) ------------------------------------
# WAIC is calculated from the full model object, which is too big for GitHub.
# (waic.msDS <- waicAbund(out.msDS, by.sp = TRUE))
# (waic.lfMsDS <- waicAbund(out.lfMsDS, by.sp = TRUE))
# (waic.sfMsDS <- waicAbund(out.sfMsDS, by.sp = TRUE))
# save(waic.msDS, waic.lfMsDS, waic.sfMsDS, file = 'results/neon-waic-results.rda')
load('results/neon-waic-results.rda')
sum(waic.msDS[, 3])
sum(waic.lfMsDS[, 3])
sum(waic.sfMsDS[, 3])

# Difference in spatial and nonspatial model without latent factors
sum(waic.sfMsDS[, 3]) - sum(waic.msDS[, 3])
# Difference in spatial and latent factor models
sum(waic.sfMsDS[, 3]) - sum(waic.lfMsDS[, 3])

# Generate maps of species-specific density -------------------------------
# Load prediction results
load("results/neon-pred-sfMsDS-results.rda")
# Load shape file for Disney Wildernees Preserve
neon <- st_read('data/neon-locations/', "terrestrialSamplingBoundaries")
dsny <- neon %>%
  filter(siteID == 'DSNY')
dsny <- dsny %>%
  st_transform(st_crs("+proj=utm +zone=17 +units=m +datum=NAD83"))

# Carolina Wren -----------------------------------------------------------
curr.sp <- which(sp.names == 'CARW')
plot.df <- data.frame(Easting = neonPredData$easting,
                      Northing = neonPredData$northing,
                      forest = neonPredData$forest,
                      grass = neonPredData$grass,
                      mu.0.mean = mu.0.means[curr.sp, ], 
                      mu.0.ci.width = mu.0.quants[3, curr.sp, ] - mu.0.quants[1, curr.sp, ])

coords.sf <- st_as_sf(x = plot.df,
                      coords = c('Easting', 'Northing'),
                      crs = "+proj=utm +zone=17 +units=m +datum=NAD83")
coords.stars <- st_as_stars(plot.df)

CARW.plot <- ggplot(data = dsny) +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.mean), interpolate = TRUE) +
  geom_sf(alpha = 0, col = 'NA') +
  scale_fill_gradientn(colors = plasma(10), na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = 'Individuals\nper ha\u00b2', x = 'Longitude', y = 'Latitude', 
       title = '(a) CARW Mean') +
  annotation_scale(text_cex=1, bar_cols = 'black', text_family = 'LM Roman 10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14), 
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 12),
        # legend.position = c(0.87, 0.85), 
        legend.background = element_rect(fill = NA))
CARW.ci.plot <- ggplot(data = dsny) +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.ci.width), 
             interpolate = TRUE) +
  geom_sf(alpha = 0, col = 'NA') +
  scale_fill_gradientn(colors = plasma(10), na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = 'Individuals\nper ha\u00b2', x = 'Longitude', y = 'Latitude', 
       title = '(c) CARW 95% CI Width') +
  annotation_scale(text_cex=1, bar_cols = 'black', text_family = 'LM Roman 10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14), 
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 12),
        # legend.position = c(0.87, 0.85), 
        legend.background = element_rect(fill = NA))
# Eastern Meadowlark ------------------------------------------------------
curr.sp <- which(sp.names == 'EAME')
plot.df <- data.frame(Easting = neonPredData$easting, 
                      Northing = neonPredData$northing,
                      mu.0.mean = mu.0.means[curr.sp, ],
                      mu.0.ci.width = mu.0.quants[3, curr.sp, ] - mu.0.quants[1, curr.sp, ])

coords.sf <- st_as_sf(x = plot.df,
                      coords = c('Easting', 'Northing'),
                      crs = "+proj=utm +zone=17 +units=m +datum=NAD83")
coords.stars <- st_as_stars(plot.df)

EAME.plot <- ggplot(data = dsny) +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.mean), interpolate = TRUE) +
  geom_sf(alpha = 0, col = 'NA') +
  scale_fill_gradientn(colors = plasma(10), na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = 'Individuals\nper ha\u00b2', x = 'Longitude', y = 'Latitude', 
       title = '(b) EAME Mean') +
  annotation_scale(text_cex=1, bar_cols = 'black', text_family = 'LM Roman 10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14), 
        legend.title = element_text(size = 12),
        text = element_text(family = 'LM Roman 10'),
        # legend.position = c(0.87, 0.85), 
        legend.background = element_rect(fill = NA))
EAME.ci.plot <- ggplot(data = dsny) +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.ci.width), 
             interpolate = TRUE) +
  geom_sf(alpha = 0, col = 'NA') +
  scale_fill_gradientn(colors = plasma(10), na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = 'Individuals\nper ha\u00b2', x = 'Longitude', y = 'Latitude', 
       title = '(d) EAME 95% CI Width') +
  annotation_scale(text_cex=1, bar_cols = 'black', text_family = 'LM Roman 10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14), 
        legend.title = element_text(size = 12),
        text = element_text(family = 'LM Roman 10'),
        # legend.position = c(0.87, 0.85), 
        legend.background = element_rect(fill = NA))

plot_grid(CARW.plot, EAME.plot, CARW.ci.plot, EAME.ci.plot, nrow = 2, ncol = 2)
# Figure S1 in the manuscript
ggsave(file = 'figures/Figure-S1.png', device = 'png', width = 9, height = 10, 
       bg = 'white')

# Generate plot of forest cover effect ------------------------------------
# Species-specific effects
beta.samples <- beta.samples.sfMsDS
beta.for.samples <- MCMCchains(beta.samples, param = 'forest', exact = FALSE)
# Community-level effects
beta.comm.samples <- beta.comm.samples.sfMsDS
beta.for.comm.samples <- MCMCchains(beta.comm.samples, param = 'forest', exact = FALSE)
# Species codes, plus the community code (COMM)
sp.names <- dimnames(neonDWP$y)[[1]]
sp.names.factor <- factor(sp.names, levels = sp.names)
sp.codes <- as.numeric(sp.names.factor)
N <- length(sp.names)
cov.plot.df <- data.frame(for.mean = apply(beta.for.samples, 2, mean),
                          for.low = apply(beta.for.samples, 2, quantile, 0.25),
                          for.lowest = apply(beta.for.samples, 2, quantile, 0.025),
                          for.high = apply(beta.for.samples, 2, quantile, 0.75),
                          for.highest = apply(beta.for.samples, 2, quantile, 0.975),
                          sp = sp.codes)
# Rearrange and add things to get the plot to display effects in increasing
# order.
cov.plot.df <- cov.plot.df %>%
  arrange(for.mean)
cov.plot.df$sp.factor <- as.character(sp.names.factor[cov.plot.df$sp])
cov.plot.df$sort.sp <- 1:N

# Add in the community level covariate
comm.plot.df <- data.frame(for.mean = mean(beta.for.comm.samples),
                           for.low = quantile(beta.for.comm.samples, 0.25),
                           for.lowest = quantile(beta.for.comm.samples, 0.025),
                           for.high = quantile(beta.for.comm.samples, 0.75),
                           for.highest = quantile(beta.for.comm.samples, 0.975),
                           sp = N + 1,
                           sp.factor = 'COMM',
                           sort.sp = N + 1)
cov.plot.df <- rbind(cov.plot.df, comm.plot.df)
cov.plot.df$sp.factor <- factor(cov.plot.df$sp.factor, levels = c(unique(cov.plot.df$sp.factor)))

# Generate Figure 1A
for.cov.plot <- ggplot(data = cov.plot.df, aes(x = sort.sp, fill = for.mean, group = sp.factor)) +
  geom_hline(yintercept = 0, col = 'black', size = 0.75, lty = 2) +
  geom_boxplot(aes(ymin = for.lowest, lower = for.low, middle = for.mean,
		   upper = for.high, ymax = for.highest), stat = 'identity', col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 14) +
  guides(fill = "none") +
  labs(x = "Species", y = "Effect of Forest Cover", title = '(a)') +
  scale_x_continuous(breaks = 1:(N+1), labels = cov.plot.df$sp.factor) +
  theme(text = element_text(family = 'LM Roman 10'), 
        axis.text.x = element_text(angle = 45, hjust = 1))
for.cov.plot


# Generate a plot of species-specific detection probability ---------------
det.int.samples <- MCMCchains(alpha.samples.sfMsDS, param = 'Intercept', exact = FALSE) 
det.means <- apply(exp(det.int.samples), 2, mean)
det.comm.means <- mean(exp(alpha.comm.samples.sfMsDS[, 1]))
det.comm.quants <- quantile(exp(alpha.comm.samples.sfMsDS[, 1]), c(0.025, 0.975))

x.vals <- seq(0, .250, length.out = 200)
N <- length(sp.names)
n.vals <- length(x.vals)
p.plot.df <- data.frame(val = NA, 
                        x.val = rep(x.vals, N), 
                        sp = rep(sp.names, each = n.vals))
for (i in 1:N) {
  indx <- ((i - 1) * n.vals + 1):(i * n.vals)
  p.plot.df$val[indx] <- gxhn(x.vals, det.means[i])
}

comm.plot.df <- data.frame(mean = gxhn(x.vals, det.comm.means), 
                           x.val = x.vals,
                           low = gxhn(x.vals, det.comm.quants[1]),
                           high = gxhn(x.vals, det.comm.quants[2]))
# Generate Figure 1B
det.plot <- ggplot(data = comm.plot.df) + 
  geom_ribbon(aes(x = x.val, ymin = low, ymax = high), fill = 'grey', 
              alpha = 0.5) +
  geom_line(data = p.plot.df, aes(x = x.val, y = val, col = sp), linewidth = 0.75, lty = 1) + 
  theme_bw(base_size = 14) +
  geom_line(aes(x = x.val, y = mean), col = 'black', linewidth = 0.75) + 
  theme(text = element_text(family = 'LM Roman 10')) +
  labs(x = 'Distance (km)', y = 'Detection Probability', col = 'Species', title = '(b)')
# Plot and save Figure 1
plot_grid(for.cov.plot, det.plot, ncol = 1)
ggsave(file = 'figures/Figure-1.png', device = 'png', width = 9, height = 10)
