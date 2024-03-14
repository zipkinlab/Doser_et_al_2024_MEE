# combine-chains.R: script to combine chains from three separate model 
#                   runs into one model output file.  
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)

# NOTE: this script will not run as the full model results are too large
#       for GitHub. Result files can be produced by running the "main-*" files, 
#       which will then allow you to run this script (changing the directories below)

# Combine objects from three chains into one ------------------------------
combine.chains <- function(out.1, out.2, out.3) {
  out <- list()
  beta.samples <- mcmc.list(out.1$beta.samples, 
                            out.2$beta.samples,
                            out.3$beta.samples)
  tau.sq.samples <- mcmc.list(out.1$tau.sq.samples, 
                              out.2$tau.sq.samples,
                              out.3$tau.sq.samples)
  out$rhat$beta <- gelman.diag(beta.samples)$psrf[, 2]
  out$rhat$tau.sq <- gelman.diag(tau.sq.samples)$psrf[, 2]
  if (class(out.1) == 'spAbund') {
    theta.samples <- mcmc.list(out.1$theta.samples, 
                               out.2$theta.samples,
                               out.3$theta.samples)
    out$rhat$theta <- gelman.diag(theta.samples)$psrf[, 2]
    out$theta.samples <- rbind(out.1$theta.samples, out.2$theta.samples, out.3$theta.samples)
    out$coords <- out.1$coords
    out$w.samples <- rbind(out.1$w.samples, out.2$w.samples, out.3$w.samples)
    out$ESS$theta <- effectiveSize(out$theta.samples)
    out$n.neighbors <- out.1$n.neighbors
    out$cov.model.indx <- out.1$cov.model.indx
    out$type <- out.1$type
  }
  if (out.1$muRE) {
    out$sigma.sq.mu.samples <- rbind(out.1$sigma.sq.mu.samples, 
                                     out.2$sigma.sq.mu.samples, 
                                     out.3$sigma.sq.mu.samples)
    out$beta.star.samples <- rbind(out.1$beta.star.samples, 
                                   out.2$beta.star.samples, 
                                   out.3$beta.star.samples)
    sigma.sq.mu.samples <- mcmc.list(out.1$sigma.sq.mu.samples, 
                                     out.2$sigma.sq.mu.samples,
                                     out.3$sigma.sq.mu.samples)
    out$rhat$sigma.sq.mu <- gelman.diag(sigma.sq.mu.samples)$psrf[, 2]
    out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
    
  }
  out$beta.samples <- rbind(out.1$beta.samples, out.2$beta.samples, out.3$beta.samples)
  out$tau.sq.samples <- rbind(out.1$tau.sq.samples, out.2$tau.sq.samples, out.3$tau.sq.samples)
  out$y.rep.samples <- rbind(out.1$y.rep.samples, out.2$y.rep.samples, out.3$y.rep.samples)
  out$like.samples <- rbind(out.1$like.samples, out.2$like.samples, out.3$like.samples)
  out$mu.samples <- rbind(out.1$mu.samples, out.2$mu.samples, out.3$mu.samples)
  out$X <- out.1$X
  out$X.re <- out.1$X.re
  out$y <- out.1$y
  out$re.level.names <- out.1$re.level.names
  out$ESS$beta <- effectiveSize(out$beta.samples)
  out$ESS$tau.sq <- effectiveSize(out$tau.sq.samples)
  out$call <- out.1$call
  out$n.samples <- out.1$n.samples
  out$n.post <- out.1$n.post
  out$n.thin <- out.1$n.thin
  out$n.burn <- out.1$n.burn
  out$n.chains <- 3
  out$re.cols <- out.1$re.cols
  out$dist <- out.1$dist
  out$muRE <- out.1$muRE
  out$run.time <- out.1$run.time
  class(out) <- class(out.1)
  out
}

# Non-spatial model -------------------------------------------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-nonspatial-chain-1.rda")
out.non.sp.1 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-nonspatial-chain-2.rda")
out.non.sp.2 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-nonspatial-chain-3.rda")
out.non.sp.3 <- out

out.non.sp <- combine.chains(out.non.sp.1, out.non.sp.2, out.non.sp.3)
save(out.non.sp, file = '/mnt/disk4/jeff/DFKZ23/results/fia-nonspatial-full.rda')
rm(out.non.sp, out.non.sp.1, out.non.sp.3, out)
gc()

# Non-spatial model + eco slopes ------------------------------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-ecoregion-chain-1.rda")
out.eco.1 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-ecoregion-chain-2.rda")
out.eco.2 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-ecoregion-chain-3.rda")
out.eco.3 <- out

out.eco <- combine.chains(out.eco.1, out.eco.2, out.eco.3)
save(out.eco, file = '/mnt/disk4/jeff/DFKZ23/results/fia-ecoregion-full.rda')
rm(out.eco, out.eco.1, out.eco.3, out)
gc()

# Spatial model -----------------------------------------------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-chain-1.rda")
out.sp.1 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-chain-2.rda")
out.sp.2 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-chain-3.rda")
out.sp.3 <- out

out.sp <- combine.chains(out.sp.1, out.sp.2, out.sp.3)
save(out.sp, file = '/mnt/disk4/jeff/DFKZ23/results/fia-sp-full.rda')
rm(out.sp, out.sp.1, out.sp.3, out)
gc()

# Spatial model + eco slopes ----------------------------------------------
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-chain-1.rda")
out.sp.eco.1 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-chain-2.rda")
out.sp.eco.2 <- out
load("/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-chain-3.rda")
out.sp.eco.3 <- out

out.sp.eco <- combine.chains(out.sp.eco.1, out.sp.eco.2, out.sp.eco.3)
save(out.sp.eco, file = '/mnt/disk4/jeff/DFKZ23/results/fia-sp-ecoregion-full.rda')
rm(out.sp.eco, out.sp.eco.1, out.sp.eco.3, out)
gc()
