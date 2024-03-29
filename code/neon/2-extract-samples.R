# 2-extract-samples.R: this script extracts MCMC samples and Rhat values from 
#                      the full multi-species model objects, which are too large
#                      to put onto GitHub.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)

load("results/neon-msDS-results.rda")
out.msDS <- out
beta.samples.msDS <- out.msDS$beta.samples
beta.comm.samples.msDS <- out.msDS$beta.comm.samples
tau.sq.beta.samples.msDS <- out.msDS$tau.sq.beta.samples
alpha.samples.msDS <- out.msDS$alpha.samples
alpha.comm.samples.msDS <- out.msDS$alpha.comm.samples
tau.sq.alpha.samples.msDS <- out.msDS$tau.sq.alpha.samples
rhat.msDS <- out.msDS$rhat
save(beta.samples.msDS, beta.comm.samples.msDS, tau.sq.beta.samples.msDS, 
     alpha.samples.msDS, alpha.comm.samples.msDS, tau.sq.alpha.samples.msDS, 
     rhat.msDS, file = 'results/neon-msDS-small-results.rda')

load("results/neon-lfMsDS-results.rda")
out.lfMsDS <- out
beta.samples.lfMsDS <- out.lfMsDS$beta.samples
beta.comm.samples.lfMsDS <- out.lfMsDS$beta.comm.samples
tau.sq.beta.samples.lfMsDS <- out.lfMsDS$tau.sq.beta.samples
alpha.samples.lfMsDS <- out.lfMsDS$alpha.samples
alpha.comm.samples.lfMsDS <- out.lfMsDS$alpha.comm.samples
tau.sq.alpha.samples.lfMsDS <- out.lfMsDS$tau.sq.alpha.samples
lambda.samples.lfMsDS <- out.lfMsDS$lambda.samples
rhat.lfMsDS <- out.lfMsDS$rhat
save(beta.samples.lfMsDS, beta.comm.samples.lfMsDS, tau.sq.beta.samples.lfMsDS, 
     alpha.samples.lfMsDS, alpha.comm.samples.lfMsDS, tau.sq.alpha.samples.lfMsDS, 
     lambda.samples.lfMsDS, rhat.lfMsDS, file = 'results/neon-lfMsDS-small-results.rda')

load("results/neon-sfMsDS-results.rda")
out.sfMsDS <- out
beta.samples.sfMsDS <- out.sfMsDS$beta.samples
beta.comm.samples.sfMsDS <- out.sfMsDS$beta.comm.samples
tau.sq.beta.samples.sfMsDS <- out.sfMsDS$tau.sq.beta.samples
alpha.samples.sfMsDS <- out.sfMsDS$alpha.samples
alpha.comm.samples.sfMsDS <- out.sfMsDS$alpha.comm.samples
tau.sq.alpha.samples.sfMsDS <- out.sfMsDS$tau.sq.alpha.samples
lambda.samples.sfMsDS <- out.sfMsDS$lambda.samples
rhat.sfMsDS <- out.sfMsDS$rhat
save(beta.samples.sfMsDS, beta.comm.samples.sfMsDS, tau.sq.beta.samples.sfMsDS, 
     alpha.samples.sfMsDS, alpha.comm.samples.sfMsDS, tau.sq.alpha.samples.sfMsDS, 
     lambda.samples.sfMsDS, rhat.sfMsDS, file = 'results/neon-sfMsDS-small-results.rda')


