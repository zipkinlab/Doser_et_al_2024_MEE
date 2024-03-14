# pred-main.R: code to predict forest AGB across the continental US using 
#              a spatial linear mixed model with a random slope of tree
#              canopy cover. This script will not run as the model object
#              is too large for GitHub.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# NOTE: this script will not run as the full model results are too large
#       for GitHub. Result files can be produced by running the "main-*" files, 
#       which will then allow you to run this script (changing the directories below)

# Read in model output object ---------------------------------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-full.rda")

# Load the data -----------------------------------------------------------
load("data/fia-data.rda")

# Read in prediction values -----------------------------------------------
load("data/fia-pred-data.rda")
# Scale prediction covariates by values used to fit model
J.0 <- nrow(pred.covs)
# The plus one is for the ecoregion random effect
X.0 <- matrix(1, J.0, ncol(out.sp.eco$X) + 1)
colnames(X.0) <- c(colnames(out.sp.eco$X), 'ecoregion')
elev.mean <- mean(data.list$covs$elev)
elev.sd <- sd(data.list$covs$elev)
tmax.mean <- mean(data.list$covs$tmax)
tmax.sd <- sd(data.list$covs$tmax)
tcc.mean <- mean(data.list$covs$tcc)
tcc.sd <- sd(data.list$covs$tcc)
X.0[, 'scale(elev)'] <- (pred.covs$elev.elevation - elev.mean) / elev.sd
X.0[, 'scale(tmax)'] <- (pred.covs$tmax - tmax.mean) / tmax.sd
X.0[, 'I(scale(tmax)^2)'] <- X.0[, 'scale(tmax)']^2 
X.0[, 'scale(tcc)'] <- (pred.covs$tcc - tcc.mean) / tcc.sd
X.0[, 'ecoregion'] <- pred.covs$ecoregion

# Predict piece by piece across the continental US. -----------------------
# Split up the data set.
vals <- split(1:J.0, ceiling(seq_along(1:J.0) / 10000))
y.0.quants <- matrix(NA, 3, J.0)
w.quants <- matrix(NA, 3, J.0)
for (j in 1:length(vals)) {
  print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[j]]
  out.pred <- predict(out.sp.eco, X.0[curr.indx, ], coords.0[curr.indx, ], n.omp.threads = 10,
		      verbose = FALSE)
  y.0.quants[, curr.indx] <- apply(out.pred$y.0.samples^2, 2, quantile, c(0.025, 0.5, 0.975))
  w.quants[, curr.indx] <- apply(out.pred$w.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
}

# Save results ------------------------------------------------------------
save(y.0.quants, w.quants, file = '/mnt/disk4/jeff/DFKZ23/results/fia-pred-results.rda')

