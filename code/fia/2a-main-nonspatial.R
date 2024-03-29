# 2a-main-nonspatial.R: script to fit a linear mixed model for the FIA biomass data
#                       where the random intercept is for ecoregion.
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
# Initial values come from an initial model run
inits <- list(beta = c(6.13, -0.74, -1.16, -0.18, 0.49),
              sigma.sq.mu = c(1.614),
              tau.sq = 5.42)
priors <- list(tau.sq.ig = c(2, 1), 
               beta.normal = list(mean = 0, var = 1000))

# Use square root of biomass to better approximate normal distribution and 
# ensure positive support on back-transformation.
data.list$y <- sqrt(data.list$y)

# Fit the model -----------------------------------------------------------
n.batch <- 10000 
batch.length <- 25
n.burn <- 190000
n.thin <- 20 
n.chains <- 1

out <- abund(formula = ~ scale(elev) + scale(tmax) + I(scale(tmax)^2) + 
                         scale(tcc) + (1 | ecoregion),
             data = data.list, priors = priors, inits = inits,
             family = 'Gaussian',
             n.batch = n.batch, batch.length = batch.length,
             n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin,
             n.chains = n.chains, n.report = 1, n.omp.threads = 1)

# Save to hard drive ------------------------------------------------------
save(out, file = paste0('/mnt/disk4/jeff/DFKZ23/results/fia-nonspatial-chain-', 
                        chain, '.rda'))
