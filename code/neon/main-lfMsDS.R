# main-lfMsDS.R: fit a multi-species distance sampling model with species
#                correlations to estimate density of 16 bird species 
#                in the Disney Wilderness Preserve. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# Set priors --------------------------------------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			               var = 10),
                   kappa.unif = list(0, 100),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

# Set initial values ------------------------------------------------------
inits.list <- list(kappa = 1, alpha.comm = 0, beta.comm = 0, alpha = 0)

# Set tuning values -------------------------------------------------------
tuning <- list(beta = 0.1, alpha = 1, beta.star = 0.3, alpha.star = 0.1,
               kappa = 0.8, lambda = 1, w = 1)

# Run the model -----------------------------------------------------------
n.batch <- 4000
batch.length <- 25
n.burn <- 50000
n.thin <- 50
n.chains <- 3

out <- lfMsDS(abund.formula = ~ scale(forest) + scale(grass),
	      det.formula = ~ scale(wind),
	      data = neonDWP,
	      n.batch = n.batch,
	      batch.length = batch.length,
	      inits = inits.list,
	      family = 'Poisson',
	      det.func = 'halfnormal',
	      transect = 'point',
	      tuning = tuning,
	      n.factors = 2,
	      priors = prior.list,
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = TRUE,
	      n.report = 10,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)

# Save results to hard drive ----------------------------------------------
save(out, file = "results/neon-lfMsDS-results.rda")
