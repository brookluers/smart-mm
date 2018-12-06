#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
cat("Running 'correct random effects', intercept-only generative, evenly spaced time points simulation with ")

library(parallel)

N <- as.numeric(args[1])
nsim <- as.numeric(args[2])
myseed <- as.numeric(args[3])
mycores <- args[4]
doEffsize <- args[5]

set.seed(myseed)
cat("\nN = "); cat(N)
cat("\nnsim = "); cat(nsim)
cat("\neffect size = "); cat(doEffsize)
if (is.na(doEffsize)) {
  doEffsize <- 'all'
}
cat("\nseed = "); cat(myseed)
cat("\ndesired cores = "); cat(mycores);
cat("\ndetectCores() = "); cat(detectCores())
mycores <- max(c(mycores, detectCores() - 1))
options(cores=mycores)
options(mc.cores=mycores)
cat("\n--using "); cat(mycores); cat(" cores")
cat("\noption('mc.cores',2) = "); cat(getOption("mc.cores", 2))

# Set up generative parameters so that the true marginal random effects structure
# is the same as the one we use for estimation.
# And control the effect size.

## Make end-of-study(time) - knot = 1

## When alhpa5=alpha6=0 and psi(a1)=0, then the
## marginal variances do not differ by regime, and the marginal residual has the 
## random effects structure we use for estimation.
##   In this case, ||Vhat - V|| should be very small
source("gen.R")
source("fit.R")
source("fitplugin.R")
# alpha <- c(1, 0.09, 0.1, 0.15, 0.08, 0, 0)
a1s <- c(1,-1)
a2s <- c(1,-1)
alpha_small <- c(1, 0.1, 0.027, 0.1, 0.08, 0, 0)
alpha_med <- c(1, 0.1, 0.125, 0.1, 0.08, 0, 0)
alpha_large <- c(1, 0.1, 0.225, 0.1, 0.08, 0, 0)
psi <- c('1' = 0, '-1' = 0)
theta <- -0.08
tvec <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
knot <- tvec[5]
sigma <- 1
cutoff <- 1.1
ff_Zgen <- Y ~ 1 # INTERCEPTS ONLY GENERATIVE MODEL 
G <- matrix(0.8)

covfunc_epsilon <- NULL

source("calc-simparm.R", print.eval=TRUE)
small_effsizes
med_effsizes
large_effsizes
cat("\nTrue variance-covariance for one individual, identical across regimes:\n")
(Vi <- Zi %*% G %*% t(Zi) + sigma^2 * diag(length(tvec)))
cat("\n")
cat("...as a correlation matrix: \n")
cov2cor(Vi)

if (doEffsize == 'small') {
  cat('\n\nsaving results for small effect size\n')
  res_small <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_small, Vtruelist$small)))
  save(alpha_small, small_effsizes, truecoefs_small, res_small, simparm,
       file= paste("sim1-smalleffect-N",N,"-nsim",nsim,".RData", sep=''))
  rm(res_small)
} else if (doEffsize == 'med') {
  cat('\n\nsaving results for medium effect size\n')
  res_med <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_med, Vtruelist$med)))
  save(alpha_med, med_effsizes, truecoefs_med, res_med, simparm,
       file= paste("sim1-medeffect-N",N,"-nsim",nsim,".RData", sep=''))
  
} else if (doEffsize == 'large') {
  cat('\n\nsaving results for large effect size\n')
  res_large <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_large, Vtruelist$large)))
  save(alpha_large, large_effsizes, truecoefs_large, res_large, simparm,
       file= paste("sim1-largeeffect-N",N,"-nsim",nsim,".RData", sep=''))
  
} else {
  cat('\nsaving results for small, medium, and large effect sizes\n')
  res_small <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_small, Vtruelist$small)))
  save(alpha_small, small_effsizes, truecoefs_small, res_small, simparm,
       file= paste("sim1b-smalleffect-N",N,"-nsim",nsim,".RData", sep=''))
  rm(res_small)
  res_med <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_med, Vtruelist$med)))
  save(alpha_med, med_effsizes, truecoefs_med, res_med, simparm,
       file= paste("sim1b-medeffect-N",N,"-nsim",nsim,".RData", sep=''))
  rm(res_med)
  res_large <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_large, Vtruelist$large)))
  save(alpha_large, large_effsizes, truecoefs_large, res_large, simparm,
       file= paste("sim1b-largeeffect-N",N,"-nsim",nsim,".RData", sep=''))
}
