# main-sp-ecoregion.R: script to fit a spatial linear mixed model with a random
#                      slope of tree canopy cover by ecoregion. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# Get chain number from command line run ----------------------------------
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# Alternatively, if not running the script from the command line:
# chain <- 1
# Or, can use the n.chains function in spAbundance for sequential runs of chains
if(length(chain) == 0) base::stop('Need to tell spAbundance the chain number')

# Load the data -----------------------------------------------------------
load("data/fia-data.rda")

# Set priors and initial values -------------------------------------------
low.dist <- 0.5
high.dist <- 1000

priors <- list(phi.unif = c(3 / high.dist, 3 / low.dist), 
               sigma.sq.ig = c(2, 1),
               tau.sq.ig = c(2, 1), 
               beta.normal = list(mean = 0, var = 1000))
# Load initial values for spatial random effects from previous model run
load("data/fia-sp-eco-inits.rda")
inits <- list(beta = c(5.87, -0.65, -1.21, -0.11, 0.55),
	      sigma.sq.mu = 0.272,
	      sigma.sq = 2.41,
	      phi = 0.0089,
	      tau.sq = 4.433,
              w = w.means)

tuning.list <- list(phi = 0.02)

# Use square root of biomass to better approximate normal distribution and 
# ensure positive support on back-transformation.
data.list$y <- sqrt(data.list$y)

# Run the model -----------------------------------------------------------
n.batch <- 10000 
batch.length <- 25
n.burn <- 190000
n.thin <- 20 
n.chains <- 1

out <- spAbund(formula = ~ scale(elev) + scale(tmax) + I(scale(tmax)^2) + 
	                   scale(tcc) + (scale(tcc) | ecoregion),
		  data = data.list, priors = priors, inits = inits,
		  tuning = tuning.list,
	          n.neighbors = 5, cov.model = 'exponential', NNGP = TRUE,
		  family = 'Gaussian',
	          n.batch = n.batch, batch.length = batch.length,
	          n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
	          n.chains = n.chains, n.report = 1, n.omp.threads = 5)

# Save to hard drive ------------------------------------------------------
save(out, file = paste0('/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-chain-', 
			chain, '.rda'))

