# main-sfMsDS.R: fit a spatial multi-species distance sampling model with species
#                correlations to estimate density of 16 bird species 
#                in the Disney Wilderness Preserve. 
rm(list = ls())
library(spAbundance)
set.seed(100)

# Set priors --------------------------------------------------------------
dist.neon <- dist(neonDWP$coords)
low.dist <- min(dist.neon)
high.dist <- max(dist.neon)
mean.dist <- mean(dist.neon)
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
		   alpha.comm.normal = list(mean = 0,
			               var = 10),
		   phi.unif = list(3 / high.dist, 3 / low.dist),
                   kappa.unif = list(0, 100),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))

# Set initial values ------------------------------------------------------
inits.list <- list(alpha.comm = c(-2.5, -0.02), tau.sq.alpha = c(0.13, 0.03),
		   beta.comm = c(-2, 0, 0), tau.sq.beta = c(1.4, 0.09, 0.05),
		   phi = c(0.0004, 0.0004))

# Set tuning values -------------------------------------------------------
tuning <- list(beta = 0.1, alpha = 0.5, beta.star = 0.3, alpha.star = 0.1,
               lambda = 0.5, w = 0.5, phi = 0.5)

# Run the model -----------------------------------------------------------
n.batch <- 4000
batch.length <- 25
n.burn <- 50000
n.thin <- 50
n.chains <- 3

out <- sfMsDS(abund.formula = ~ scale(forest) + scale(grass),
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
	      NNGP = TRUE, 
	      n.neighbors = 15,
	      priors = prior.list,
	      accept.rate = 0.43,
	      n.omp.threads = 1,
	      verbose = TRUE,
	      n.report = 10,
	      n.burn = n.burn,
	      n.thin = n.thin,
	      n.chains = n.chains)

# Save results to hard drive ----------------------------------------------
save(out, file = "results/neon-sfMsDS-results.rda")
