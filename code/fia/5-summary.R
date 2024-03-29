# 5-summary.R: this script summarizes results from the FIA biomass case study
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)
library(tidyverse)
library(sf)
library(viridis)
library(MCMCvis)
library(stars)
library(patchwork)

# Load the data -----------------------------------------------------------
load("data/fia-data.rda")

# Load the results --------------------------------------------------------
# WAIC results
load("results/fia-waic-results.rda")
# Difference in WAIC between top model and all candidate models
waic.vals[4] - waic.vals
load("results/fia-top-model-results.rda")

# Plot of fitted values vs. the true values -------------------------------
# Not surprisingly, the model does not do great at predicting extremely high 
# biomass values.
y.true <- data.list$y
plot(y.true, y.rep.quants[2, ], pch = 19)
abline(0, 1)

# TCC effects -------------------------
beta.tcc.means <- apply(beta.star.samples + beta.samples[, 5], 2, median)

# Stuff for plots ---------------------------------------------------------
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
ecr <- st_read(dsn = "data/usgs-ecoregions/", layer = "us_eco_l3")
ecrs.albers <- ecr %>%
  st_transform(crs = my.proj)
ecrs.albers$beta.tcc <- rep(NA, nrow(ecrs.albers))
for (j in 1:nrow(ecrs.albers)) {
  ecrs.albers$beta.tcc[j] <- beta.tcc.means[as.numeric(ecrs.albers$US_L3CODE[j])]
}
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = my.proj)

# Generate Figure 2 -------------------------------------------------------
load("data/fia-pred-data.rda")
coords.pred.sf <- st_as_sf(as.data.frame(coords.0), 
                           coords = c('X', 'Y'),
                           crs = my.proj)
load("results/fia-pred-results.rda")
pred.plot.df <- data.frame(x = coords.0[, 1],
                           y = coords.0[, 2],
                           y.med = y.0.quants[2, ],
                           tcc = pred.covs$tcc,
                           y.ci.width = y.0.quants[3, ] - y.0.quants[1, ],
                           w.med = w.quants[2, ],
                           w.ci.width = w.quants[3, ] - w.quants[2, ])
pred.stars.df <- st_as_stars(pred.plot.df, dims = c('x', 'y'))

coords.sf <- st_as_sf(as.data.frame(data.list$coords), 
                      coords = c('x', 'y'),
                      crs = my.proj)
# Figure 2A ----------------------------
points.plot <- ggplot(coords.sf) +
  geom_sf(size = 0.005, col = 'black') +
  geom_sf(data = usa, fill = NA, color=alpha("gray", 0.75), lwd = 0.6) +
  theme_bw(base_size = 18) +
  labs(title = '(a) Observed locations') +
  theme(text = element_text(family="LM Roman 10"),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
	plot.title = element_text(size = 14))

# Figure 2B ---------------------------
beta.tcc.plot <- ggplot(ecrs.albers) + 
  geom_sf(aes(fill = beta.tcc)) + 
  theme_bw(base_size = 14) + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
    	               na.value = NA) + 
  labs(title = '(b) TCC median effect', fill = '') + 
  theme(legend.position = c(0.92, 0.28), 
        plot.title = element_text(size = 14), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'))

# Figure 2C ---------------------------
y.med.plot <- ggplot() +
  geom_stars(data = pred.stars.df, aes(x = x, y = y, fill = y.med),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = plasma(10), na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = "", title = '(c) Median biomass') +
  theme(legend.position = c(0.92, 0.28), 
        plot.title = element_text(size = 14), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'))

# Figure 2D ---------------------------
y.ci.width.plot <- ggplot() +
  geom_stars(data = pred.stars.df, aes(x = x, y = y, fill = y.ci.width),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = plasma(10), na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = "", title = '(d) 95% CI for biomass') +
  theme(legend.position = c(0.92, 0.28), 
        plot.title = element_text(size = 14), 
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'))

# Figure 2
my.plot <- points.plot + beta.tcc.plot + y.med.plot + y.ci.width.plot
ggsave(my.plot, file = 'figures/Figure-2.png', width = 10, height = 7, units = 'in', 
       bg = 'white')
