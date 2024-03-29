# 2b-main-ecoregion.R: script to fit a linear mixed model with a random
#                      slope of tree canopy cover by ecoregion and a random
#                      intercept of ecoregion.. 
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
inits <- list(beta = c(6.24, -0.35, -0.83, -0.15, 0.67),
              sigma.sq.mu = c(0.41, 1.72),
              tau.sq = c(5.22))
priors <- list(tau.sq.ig = c(2, 1), 
               beta.normal = list(mean = 0, var = 1000))

# Use square root of biomass to better approximate normal distribution and 
# ensure positive support on back-transformation.
data.list$y <- sqrt(data.list$y)

# Run the model -----------------------------------------------------------
n.batch <- 10000 
batch.length <- 25
n.burn <- 190000
n.thin <- 20 
n.chains <- 1

out <- abund(formula = ~ scale(elev) + scale(tmax) + I(scale(tmax)^2) + 
                         scale(tcc) + (scale(tcc) | ecoregion) + (1 | ecoregion),
             data = data.list, priors = priors, inits = inits,
             family = 'Gaussian',
             n.batch = n.batch, batch.length = batch.length,
             n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
             n.chains = n.chains, n.report = 1, n.omp.threads = 1)

# Save the model object ---------------------------------------------------
save(out, file = paste0('/mnt/disk4/jeff/DFKZ23/results/fia-ecoregion-chain-', 
                        chain, '.rda'))
