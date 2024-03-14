# predict-sfMsDS.R: this script predict density for the sixteen bird species
#                   across the Disney Wilderness Preserve using results from 
#                   a spatially-explicit hierarchical distance sampling model.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# Read in data and fitted model -------------------------------------------
load("results/neon-sfMsDS-results.rda")
# Read in prediction data -------------------------------------------------
grass.pred <- (neonPredData[, 'grass'] - mean(neonDWP$covs$grass)) / sd(neonDWP$covs$grass)
forest.pred <- (neonPredData[, 'forest'] - mean(neonDWP$covs$forest)) / sd(neonDWP$covs$forest)
X.0 <- cbind(1, forest.pred, grass.pred)
colnames(X.0) <- c('(Intercept)', 'scale(forest)', 'scale(grass)')
coords.0 <- as.matrix(neonPredData[, c('easting', 'northing')])
pred.out <- predict(out, X.0, coords.0, n.omp.threads = 1)
mu.0.means <- apply(pred.out$mu.0.samples, c(2, 3), mean)
mu.0.quants <- apply(pred.out$mu.0.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
save(mu.0.means, mu.0.quants, file = 'results/neon-pred-sfMsDS-results.rda')
