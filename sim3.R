#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
cat("Running 'simulation 3' with ")

library(parallel)

N <- as.numeric(args[1])
nsim <- as.numeric(args[2])
myseed <- as.numeric(args[3])
mycores <- args[4]
doEffsize <- args[5]
missprob <- args[6]
set.seed(myseed)
cat("\nN = "); cat(N)
cat("\nnsim = "); cat(nsim)
cat("\neffect size = "); cat(doEffsize)
if (is.na(doEffsize)) {
  doEffsize <- 'all'
}
if (is.na(missprob)){
  cat("No missingness probability specified, default to 0.2\n")
  missprob <- 0.2
} else {
  missprob <- as.numeric(missprob)
}

cat("\nseed = "); cat(myseed)
cat("\ndesired cores = "); cat(mycores);
cat("\ndetectCores() = "); cat(detectCores())
mycores <- max(c(mycores, detectCores() - 1))
options(cores=mycores)
options(mc.cores=mycores)
cat("\n--using "); cat(mycores); cat(" cores")
cat("\noption('mc.cores',2) = "); cat(getOption("mc.cores", 2))
cat("\nMissingness probability = "); cat(missprob); cat("\n")
cat("class(missprob) = "); cat(class(missprob))
source("gen.R")
source("fit.R")
source("fitplugin.R")
a1s <- c(1,-1)
a2s <- c(1,-1)
alpha_small <- c(1, 0.1, 0.057, 0.15, 0.08, 0.25, -0.27)
alpha_med <- c(1, 0.1, 0.285, 0.15, 0.08, 0.25, -0.27)
alpha_large <- c(1, 0.1, 0.512, 0.15, 0.08, 0.25, -0.27)
psi <- c('1'= 0.1, '-1' = -0.1)
theta <- -0.2
tvec <- c(0, 0.5, 1.5, 2, 2.25, 2.5, 3)
knot <- tvec[4]
sigma <- 1
cutoff <- 1.1
ff_Zgen <- Y ~ 1 + time
G <- matrix(c(0.8, -0.2, -0.2, 1),
            nrow=2,byrow=T)


missFunc <- function(dd){
  ret <- dd[sample(c(TRUE,FALSE),
            size=nrow(dd), 
            replace = TRUE, prob = c(1-missprob, missprob)),]
  nii <- table(ret[,'id'])
  names(nii)[which(nii<2)]
  return(
    ret[!(ret[,'id'] %in% as.numeric(names(nii)[which(nii<2)])),] # throw out ids where ni < 2
  )
}
covfunc_epsilon <- NULL

source("calc-simparm.R", print.eval=TRUE)

if (doEffsize == 'small') {
  cat('\n\nsaving results for small effect size\n')
  res_small <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_small, Vtruelist$small, missFunc = missFunc)))
  save(alpha_small, small_effsizes, truecoefs_small, res_small, simparm,
       file= paste("sim3-smalleffect-N",N,"-nsim",nsim,".RData", sep=''))
  rm(res_small)
} else if (doEffsize == 'med') { 
  cat('\n\nsaving results for medium effect size\n')
  res_med <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_med, Vtruelist$med, missFunc = missFunc)))
  save(alpha_med, med_effsizes, truecoefs_med, res_med, simparm,
       file= paste("sim3-medeffect-N",N,"-nsim",nsim,".RData", sep=''))
  
} else if (doEffsize == 'large') {
  cat('\n\nsaving results for large effect size\n')
  res_large <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_large, Vtruelist$large, missFunc = missFunc)))
  save(alpha_large, large_effsizes, truecoefs_large, res_large, simparm,
       file= paste("sim3-largeeffect-N",N,"-nsim",nsim,".RData", sep=''))
} else {
  cat('\nsaving results for small, medium, and large effect sizes\n')
  res_small <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_small, Vtruelist$small, missFunc = missFunc)))
  save(alpha_small, small_effsizes, truecoefs_small, res_small, simparm,
       file= paste("sim3-smalleffect-N",N,"-nsim",nsim,".RData", sep=''))
  rm(res_small)
  res_med <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_med, Vtruelist$med, missFunc = missFunc)))
  save(alpha_med, med_effsizes, truecoefs_med, res_med, simparm,
       file= paste("sim3-medeffect-N",N,"-nsim",nsim,".RData", sep=''))
  rm(res_med)
  res_large <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_large, Vtruelist$large, missFunc = missFunc)))
  save(alpha_large, large_effsizes, truecoefs_large, res_large, simparm,
       file= paste("sim3-largeeffect-N",N,"-nsim",nsim,".RData", sep=''))
}
