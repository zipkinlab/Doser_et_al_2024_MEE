# 2-predict-NMix.R: script to predict abundance across Hubbard Brook Experimental
#                   Forest for the Black-throated Blue Warbler using a non-spatial
#                   Poisson N-mixture model (the top performing model).
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(spAbundance)
# Extract HBEF data for BTBW from spAbudance ------------------------------
data.hbef <- hbefCount2015
data.hbef$y <- data.hbef$y[which(dimnames(hbefCount2015$y)[[1]] == 'BTBW'), , ]

# Read in results --------------------------------------------------------- 
# Poisson N-mixture -------------------
load('results/hbef-NMix-poisson-fit.rda')
out.nmix <- out

# Predict abundance across the region -------------------------------------
# Standardize elevation values by those used to fit the model. 
# Note the hbefElev values come from the spOccupancy R package. 
elev.pred <- (hbefElev$val - mean(data.hbef$abund.covs[, 1])) / sd(data.hbef$abund.covs[, 1])
# This is the design matrix for prediction.
X.0 <- cbind(1, elev.pred, elev.pred^2)
# Generate predictions across the region.
out.nmix.pred <- predict(out.nmix, X.0)
# Expected abundance means
mu.0.nmix.mean <- apply(out.nmix.pred$mu.0.samples, 2, mean)
# Expected abundance SDs
mu.0.nmix.sd <- apply(out.nmix.pred$mu.0.samples, 2, sd)
# Prediction coordinates
coords.0 <- as.matrix(hbefElev[, c('Easting', 'Northing')])

save(mu.0.nmix.mean, mu.0.nmix.sd, coords.0, 
     file = 'results/hbef-NMix-poisson-pred-results.rda')
