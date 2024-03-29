# 1c-main-spNMix.R: fit a spatial N-mixture model with a Poisson distribution
#                   to estimate abundance of Black-throated Blue Warblers across
#                   Hubbard Brook Experimental Forest.
# Author: Jeffrey W. Doser
rm(list = ls())
set.seed(400)
library(spAbundance)

# Extract HBEF data for BTBW from spAbudance ------------------------------
data.hbef <- hbefCount2015
data.hbef$y <- data.hbef$y[which(dimnames(hbefCount2015$y)[[1]] == 'BTBW'), , ]

# Specify priors, initial values, tuning values for spAbundance -----------
# Priors ------------------------------
dist.coords <- dist(data.hbef$coords)
low.dist <- quantile(dist.coords, 0.01)
mean.dist <- mean(dist.coords)
max.dist <- max(dist.coords)
prior.list <- list(beta.normal = list(0, 100),
                   alpha.normal = list(0, 2.72),
                   phi = c(3 / max.dist, 3 / low.dist), 
                   sigma.sq = c(2, 1))
# Starting values ---------------------
inits.list <- list(alpha = 0, beta = 0, w = rep(0, nrow(data.hbef$coords)),
                   phi = 3 / mean.dist, sigma.sq = 1,
                   N = apply(data.hbef$y, 1, max, na.rm = TRUE))
# Tuning values -----------------------
# Good starting values for the tuning parameters would be the estimated 
# standard deviation of the resulting estimates. 
tuning.list <- list(phi = 0.5, beta = 0.1, alpha = 0.1, w = 0.5, kappa = 0.2)

# Fit the model -----------------------------------------------------------
n.batch <- 5000
batch.length <- 25
n.burn <- 85000
n.thin <- 20
n.chains <- 3

out <- spNMix(abund.formula = ~ scale(elev) + I(scale(elev)^2),
              det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod),
              data = data.hbef,
              n.batch = n.batch,
              batch.length = batch.length,
              tuning = tuning.list,
              inits = inits.list,
              priors = prior.list,
              cov.model = 'exponential',
              NNGP = TRUE,
              n.neighbors = 5,
              family = 'Poisson',
              accept.rate = 0.43,
              n.omp.threads = 1,
              verbose = TRUE,
              n.report = 10,
              n.burn = n.burn,
              n.thin = n.thin,
              n.chains = n.chains)

# Save results
save(out, file = 'results/hbef-spNMix-poisson-fit.rda')

