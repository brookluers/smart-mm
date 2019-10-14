#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
source("simfuncs.R")

a1s <- c(1,-1)
a2s <- c(1,-1)
alpha_small <- c(1, 0.1, 0.067, 0.15, 0.1, 0.25, -0.4)
alpha_med <- c(1, 0.1, 0.37, 0.15, 0.1, 0.25, -0.4)
alpha_large <- c(1, 0.1, 0.68, 0.15, 0.1, 0.25, -0.4)
effsizenames <- c('small','med','large')
alphalist <- setNames(list(alpha_small, alpha_med, alpha_large),
                      effsizenames)
psi <- c('1'= 0.1, '-1' = -0.1)
theta <- -0.2
tvec <- c(0, 0.5, 1.5, 2, 2.25, 2.5, 3)
knot <- tvec[4]
sigma <- 1
cutoff <- 1.1
ff_Zgen <- Y ~ 1 + time
G <- matrix(c(1, -0.4, -0.4, 2),
            nrow=2,byrow=T)
covfunc_epsilon <- NULL

simparm <- get_simparm(args, a1s, a2s, alphalist, effsizenames,
                        psi, theta, tvec, knot, sigma, cutoff, 
                        ff_Zgen, G, covfunc_epsilon)

cat("effect sizes: \n")
print(simparm$effsizelist)

runsim(simparm, 'sim2')