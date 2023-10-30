# summary.R: summarize results from analysis of black-throated blue warbler
#            abundance using N-mixture models in spAbundance
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(stars)
library(viridis)
library(spAbundance)
library(spOccupancy)
library(cowplot)

# Extract HBEF data for BTBW from spAbudance ------------------------------
data.hbef <- hbefCount2015
data.hbef$y <- data.hbef$y[which(dimnames(hbefCount2015$y)[[1]] == 'BTBW'), , ]

# Read in results --------------------------------------------------------- 
# Poisson N-mixture -------------------
load('results/hbef-NMix-poisson-fit.rda')
out.nmix <- out
# NB N-mixture ------------------------
load('results/hbef-NMix-NB-fit.rda')
out.nmix.nb <- out
# Spatial Poisson N-mixture -----------
load('results/hbef-spNMix-poisson-fit.rda')
out.sp.nmix <- out
# Spatial NB N-mixture ----------------
load('results/hbef-spNMix-NB-fit.rda')
out.sp.nmix.nb <- out

# Model comparison --------------------------------------------------------
waic.nmix <- waicAbund(out.nmix)
waic.nmix.nb <- waicAbund(out.nmix.nb)
waic.sp.nmix <- waicAbund(out.sp.nmix)
waic.sp.nmix.nb <- waicAbund(out.sp.nmix.nb)
# Non-spatial, Poisson N-mixture is best according to WAIC. 
waic.nmix
waic.nmix.nb
waic.sp.nmix
waic.sp.nmix.nb

# Look at results from top performing model
summary(out.nmix)

# Posterior predictive checks for top performing model --------------------
ppc.nmix <- ppcAbund(out.nmix, fit.stat = 'freeman-tukey', group = 0)
ppc.nmix.sites <- ppcAbund(out.nmix, fit.stat = 'freeman-tukey', group = 1)
ppc.nmix.visit <- ppcAbund(out.nmix, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.nmix)
summary(ppc.nmix.sites)
summary(ppc.nmix.visit)

# Generate prediction map for top model -----------------------------------
# Non-spatial N-mixture with Poisson distribution for abundance
load('results/hbef-NMix-poisson-pred-results.rda')
plot.df <- data.frame(Easting = hbefElev$Easting,
		      Northing = hbefElev$Northing,
		      elev = hbefElev$val,
		      mu.0.nmix.mean = mu.0.nmix.mean,
		      mu.0.nmix.sd = mu.0.nmix.sd)

coords.sf <- st_as_sf(x = plot.df,
		      coords = c('Easting', 'Northing'),
		      crs = "+proj=utm +zone=19 +units=m +datum=NAD83")

# Read in a shapefile of the HBEF watersheds to get the outline of the forest.
hbef <- st_read('data/hbef-spatial/', 'hbef_wsheds')
# Join all watersheds together
hbef <- st_union(hbef)
coords.stars <- st_as_stars(plot.df)

# N-mixture expected abundance
nmix.mean.plot <- ggplot(data = hbef) +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.nmix.mean)) +
  geom_sf(alpha = 0, col = 'gray58') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = 'Expected\nAbundance', x = 'Longitude', y = 'Latitude', 
       title = '(A) Posterior Mean') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(family="LM Roman 10"),
        plot.title = element_text(size = 14))

# N-mixture standard deviation
nmix.sd.plot <- ggplot(data = hbef) +
  geom_stars(data = coords.stars, aes(x = Easting, y = Northing, fill = mu.0.nmix.sd)) +
  geom_sf(alpha = 0, col = 'gray58') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 14) +
  labs(fill = 'Expected\nAbundance', x = 'Longitude', y = 'Latitude', 
       title = '(B) Standard Deviation') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(family="LM Roman 10"),
        plot.title = element_text(size = 14))

# Generate Figure S2 in the manuscript
plot_grid(nmix.mean.plot, nmix.sd.plot, ncol = 1)
# Save figure to hard drive if desired.
ggsave("figures/Figure-S2.png", height = 8, width = 7, units = 'in', bg = 'white')

