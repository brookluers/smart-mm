#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
source("simfuncs.R")

# Set up generative parameters so that the true marginal random effects structure
# is the same as the one we use for estimation.
# And control the effect size.

## Make end-of-study(time) - knot = 1

## When alhpa5=alpha6=0 and psi(a1)=0, then the
## marginal variances do not differ by regime, and the marginal residual has the 
## random effects structure we use for estimation.
##   In this case, ||Vhat - V|| should be very smal

a1s <- c(1,-1)
a2s <- c(1,-1)
alpha_small <- c(1, 0.1, 0.11, 0.15, 0.08, 0, 0)
alpha_med <- c(1, 0.1, 0.35, 0.15, 0.08, 0, 0)
alpha_large <- c(1, 0.1, 0.58, 0.15, 0.08, 0, 0)
effsizenames <- c('small','med','large')
alphalist <- setNames(list(alpha_small, alpha_med, alpha_large),
                      effsizenames)
psi <- c('1'=0, '-1' = 0)
theta <- -0.08
tvec <- c(0, 0.5, 1.5, 2, 2.25, 2.5, 3)
knot <- tvec[4]
sigma <- 1
cutoff <- 1.1
ff_Zgen <- Y ~ 1 + time
G <- matrix(c(0.8, -0.2, -0.2, 1),
            nrow=2,byrow=T)
covfunc_epsilon <- NULL

simparm <- get_simparm(args, N, a1s, a2s, alphalist, effsizenames,
                        psi, theta, tvec, knot, sigma, cutoff, 
                        ff_Zgen, G, covfunc_epsilon)

runsim(simparm, 'sim1')
