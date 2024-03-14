# get-waic-and-samples.rda: this script extracts the WAIC for each of the 
#                           four candidate models. This script will not 
#                           run on GitHub, as the model objects are too 
#                           large. The objects created by the script are
#                           available on GitHub.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)

# NOTE: this script will not run as the full model results are too large
#       for GitHub. Result files can be produced by running the "main-*" files, 
#       which will then allow you to run this script (changing the directories below)

# Load model results ------------------------------------------------------
# These are not available on GitHub.
# Non-spatial model -------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-nonspatial-full.rda")
# Non-spatial model + eco slopes ------
load("/mnt/disk4/jeff/DFKZ23/results/fia-ecoregion-full.rda")
# Spatial model -----------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-full.rda")
# Spatial model + eco sloeps ----------
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-full.rda")

# Compare models with WAIC ------------------------------------------------
waic.vals <- rep(NA, 4)
names(waic.vals) <- c('nonspatial', 'eco', 'spatial', 'spatialEco')
waic.vals[1] <- waicAbund(out.non.sp)[3]
waic.vals[2] <- waicAbund(out.eco)[3]
waic.vals[3] <- waicAbund(out.sp)[3]
waic.vals[4] <- waicAbund(out.sp.eco)[3]
# Save WAIC values to include on GitHub.
save(waic.vals, file = 'results/fia-waic-results.rda')

# Save certain parameters from top-performing model -----------------------
beta.samples <- out.sp.eco$beta.samples
beta.star.samples <- out.sp.eco$beta.star.samples
sigma.sq.mu.samples <- out.sp.eco$sigma.sq.mu.samples
theta.samples <- out.sp.eco$theta.samples
tau.sq.samples <- out.sp.eco$tau.sq.samples
y.rep.quants <- apply(out.sp.eco$y.rep.samples^2, 2, quantile, c(0.025, 0.5, 0.975))
y.rep.means <- apply(out.sp.eco$y.rep.samples^2, 2, mean)
y.rep.sd <- apply(out.sp.eco$y.rep.samples^2, 2, sd)

# Save smaller object to include on GitHub.
save(beta.samples, beta.star.samples, sigma.sq.mu.samples, 
     theta.samples, tau.sq.samples, y.rep.quants, y.rep.means, 
     y.rep.sd, file = 'results/fia-top-model-results.rda')
