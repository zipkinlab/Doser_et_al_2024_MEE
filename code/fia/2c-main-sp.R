# 2c-main-sp.R: script to fit a spatial linear mixed model.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# Get chain number from command line run ----------------------------------
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# If not running the script from the command line, set the chain number manually:
# NOTE: comment this line if running from the command line.
chain <- 1
# Or, can use the n.chains function in spAbundance for sequential runs of chains
if(length(chain) == 0) base::stop('Need to tell spAbundance the chain number')

# Load the data -----------------------------------------------------------
load("data/fia-data.rda")

# Set priors and initial values -------------------------------------------
low.dist <- 0.5
high.dist <- 550

priors <- list(phi.unif = c(3 / high.dist, 3 / low.dist), 
               sigma.sq.ig = c(2, 1),
               tau.sq.ig = c(2, 1), 
               beta.normal = list(mean = 0, var = 1000))
# Initial values based on previous model run
inits <- list(beta = c(5.66, -0.79, -1.45, -0.16, 0.36),
              sigma.sq = 2.85,
              phi = 0.0093,
              tau.sq = 4.456) 
tuning.list <- list(phi = 0.03)

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
                           scale(tcc),
               data = data.list, priors = priors, inits = inits,
               tuning = tuning.list,
               n.neighbors = 5, cov.model = 'exponential', NNGP = TRUE,
               family = 'Gaussian',
               n.batch = n.batch, batch.length = batch.length,
               n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
               n.chains = n.chains, n.report = 1, n.omp.threads = 5)

# Save to hard drive ------------------------------------------------------
save(out, file = paste0('/mnt/disk4/jeff/DFKZ23/results/fia-sp-chain-', 
                        chain, '.rda'))

