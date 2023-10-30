# main-NMix-NB.R: fit a non-spatial N-mixture model with a negative binomial
#                 distribution to estimate abundance of Black-throated Blue
#                 Warblers across Hubbard Brook Experimental Forest.
# Author: Jeffrey W. Doser
rm(list = ls())
set.seed(400)
library(spAbundance)

# Extract HBEF data for BTBW from spAbudance ------------------------------
data.hbef <- hbefCount2015
data.hbef$y <- data.hbef$y[which(dimnames(hbefCount2015$y)[[1]] == 'BTBW'), , ]

# Specify priors, initial values, tuning values for spAbundance -----------
# Priors ------------------------------
prior.list <- list(beta.normal = list(0, 100),
		   alpha.normal = list(0, 2.72),
                   kappa.unif = c(0, 100))
# Starting values ---------------------
inits.list <- list(alpha = 0, beta = 0, kappa = 1, 
		   N = apply(data.hbef$y, 1, max, na.rm = TRUE))
# Tuning values -----------------------
# Good starting values for the tuning parameters would be the estimated 
# standard deviation of the resulting estimates. 
tuning.list <- list(beta = 0.1, alpha = 0.1, kappa = 0.2)

# Fit the model -----------------------------------------------------------
# Big
n.batch <- 5000
batch.length <- 25
n.burn <- 85000
n.thin <- 20
n.chains <- 3

out <- NMix(abund.formula = ~ scale(elev) + I(scale(elev)^2),
            det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod),
            data = data.hbef,
            n.batch = n.batch,
            batch.length = batch.length,
            tuning = tuning.list,
            inits = inits.list,
            priors = prior.list,
            family = 'NB',
            accept.rate = 0.43,
            n.omp.threads = 1,
            verbose = TRUE,
            n.report = 100,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains)

# Save results
save(out, file = 'results/hbef-NMix-NB-fit.rda')
