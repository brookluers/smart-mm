#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
cat("Running the 'correct random effects' simulation with ")

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
allregime <- expand.grid(a1=a1s,a2=a2s)
regimenames <- as.character(with(allregime,interaction(a1,a2)))
alpha_small <- c(1, 0.1, 0.027, 0.15, 0.08, 0, 0)
alpha_med <- c(1, 0.1, 0.128, 0.15, 0.08, 0, 0)
alpha_large <- c(1, 0.1, 0.229, 0.15, 0.08, 0, 0)
psi <- c('1'=0, '-1' = 0)
theta <- -0.08
coefnames <- c(paste('beta', 0:6,sep=''), 'eta')
tvec <- c(0, 0.5, 1.5, 2, 2.25, 2.5, 3)
knot <- tvec[4]
sigma <- 1
X <- vector('numeric', length=N)
X[1:floor(N/2)] <- 1
X[(floor(N/2)+1):N] <- -1
cutoff <- 1.1
ff_Zgen <- Y ~ 1 + time
Zi <- model.matrix(ff_Zgen, data.frame(time=tvec,Y=-99))
G <- matrix(c(0.2,-0.05,-0.05, 0.1),
            nrow=2,byrow=T)
cat("True marginal V:\n")
(Vi <- Zi %*% G %*% t(Zi) + sigma^2 * diag(length(tvec)))
mmean <- get_mmean(alpha_small, unique(X), theta, knot, tvec, a1=c(1,-1), a2=c(1,-1), G,ff_Zgen, sigma, cutoff)
cat("\nTrue regression coefficients, small effect size: \n")
(truecoefs_small <- get_truecoefs(alpha_small, theta, knot, G, ff_Zgen, sigma, cutoff))
truecoefs_med <- get_truecoefs(alpha_med, theta, knot, G, ff_Zgen, sigma, cutoff)
truecoefs_large <- get_truecoefs(alpha_large, theta, knot, G, ff_Zgen, sigma, cutoff)

cat("\nProbability of nonresponder: \n")
(pinr_a1 <- c('1' = get_pinr(alpha_small, knot, a1=1, G, ff_Zgen, sigma, cutoff),
             '-1' = get_pinr(alpha_small, knot, a1=-1, G, ff_Zgen, sigma, cutoff)) )

vardat <- expand.grid(tval=tvec,
                      a1=c(1,-1),
                      a2=c(1,-1))
vardat$mvar <- 
  mapply(tval = vardat$tval,
         a1=vardat$a1,
         a2=vardat$a2,
         FUN = function(tval, a1, a2) return(get_margvar(tval =tval, 
                                                         alpha_small, psi, knot, 
                                                         a1=a1, a2=a2, G, ff_Zgen, sigma,cutoff)))

# All regimes have the same variance since alpha5=alpha6=psi=0
#Vtrue_a1a2 <- mapply(a1=allregime$a1, a2=allregime$a2,
#                     function(a1, a2) return(get_vcovmat(a1=a1,a2=a2,tvec,alpha_small,theta,psi,knot,G,ff_Zgen,sigma,cutoff)),
#                     SIMPLIFY = FALSE)
#names(Vtrue_a1a2) <- regimenames

ff_po <- 
  Ya1a2 ~
  1 + 
  I( time * 1 * (time <= knot) + knot * 1 * (time > knot)) +
  I( time * 1 * (time <= knot) * a1 + knot * a1 * 1 * (time > knot)) +
  I( (time - knot) * 1 * (time > knot)) + 
  I( (time - knot) * a1 * 1 * (time > knot)) +
  I( (time - knot) * a2 * 1 * (time > knot)) + 
  I( (time - knot) * a1 * a2 * 1 * (time > knot)) +
  X 

ff_Z_slopes <- Y ~ 1 + time

ff_Z_intercept <- Y ~ 1

ff_lmer_slopes <-
  Y ~
  1 + 
  I( time * 1 * (time <= knot) + knot * 1 * (time > knot)) +
  I( time * 1 * (time <= knot) * A1 + knot * A1 * 1 * (time > knot)) +
  I( (time - knot) * 1 * (time > knot)) + 
  I( (time - knot) * A1 * 1 * (time > knot)) +
  I( (time - knot) * A2 * 1 * (time > knot)) + 
  I( (time - knot) * A1 * A2 * 1 * (time > knot)) +
  X + 
  (1 + time| id2rep)

ff_lmer_intercept <-
  Y ~
  1 + 
  I( time * 1 * (time <= knot) + knot * 1 * (time > knot)) +
  I( time * 1 * (time <= knot) * A1 + knot * A1 * 1 * (time > knot)) +
  I( (time - knot) * 1 * (time > knot)) + 
  I( (time - knot) * A1 * 1 * (time > knot)) +
  I( (time - knot) * A2 * 1 * (time > knot)) + 
  I( (time - knot) * A1 * A2 * 1 * (time > knot)) +
  X + 
  (1 | id2rep)

ff_fixef <- 
  Y ~
  1 + 
  I( time * 1 * (time <= knot) + knot * 1 * (time > knot)) +
  I( time * 1 * (time <= knot) * A1 + knot * A1 * 1 * (time > knot)) +
  I( (time - knot) * 1 * (time > knot)) + 
  I( (time - knot) * A1 * 1 * (time > knot)) +
  I( (time - knot) * A2 * 1 * (time > knot)) + 
  I( (time - knot) * A1 * A2 * 1 * (time > knot)) +
  X 

ni <- length(tvec)
regimes_effnames <- expand.grid(effsize=c('small','med','large'),
                                a1=a1s, a2=a2s)
Vtruelist <- 
  lapply(c('small','med','large'), function(effsize) {
    ret <-
      mapply(a1=allregime$a1, a2= allregime$a2,
             FUN = function(a1, a2){
               SigmaTrueX1 <- matrix(0, nrow=ni, ncol=ni)
               SigmaTrueX1 <- 
                 matrix(mapply(t_ix = row(SigmaTrueX1),
                               s_ix = col(SigmaTrueX1),
                               FUN = function(t_ix, s_ix) {
                                 cij <- get_margcov(tval=tvec[t_ix], sval=tvec[s_ix],
                                                    alpha=get(paste('alpha_', effsize, sep='')), 
                                                    psi=psi, knot=knot, a1=a1, a2=a2,
                                                    G=G, ff_Z=ff_Zgen, sigma=sigma,cutoff=cutoff)
                                 return(
                                   cij
                                 )
                               }),
                        nrow=ni,byrow=T)
               return(SigmaTrueX1)
             },
             SIMPLIFY = FALSE) 
    names(ret) <- regimenames
    return(ret)
  })
names(Vtruelist) <- c('small','med','large')

small_effsizes <- 
  rbind(
  get_effsizes(a11=1, a21=1, a12=-1, a22=1, alpha= alpha_small, 
              psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
  get_effsizes(a11=1, a21=1, a12=-1, a22=-1, alpha= alpha_small, 
             psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
get_effsizes(a11=1, a21=-1, a12=-1, a22=1, alpha= alpha_small, 
             psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
get_effsizes(a11=-1, a21=-1, a12=-1, a22=1, alpha= alpha_small, 
             psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff))


med_effsizes <- 
  rbind(
    get_effsizes(a11=1, a21=1, a12=-1, a22=1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
    get_effsizes(a11=1, a21=1, a12=-1, a22=-1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
    get_effsizes(a11=1, a21=-1, a12=-1, a22=1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
    get_effsizes(a11=-1, a21=-1, a12=-1, a22=1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff))



large_effsizes <- 
  rbind(
    get_effsizes(a11=1, a21=1, a12=-1, a22=1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
    get_effsizes(a11=1, a21=1, a12=-1, a22=-1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
    get_effsizes(a11=1, a21=-1, a12=-1, a22=1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff),
    get_effsizes(a11=-1, a21=-1, a12=-1, a22=1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff))

vcomplist_a1a2_true <- setNames(lapply(regimenames,function(x)return(list(Vhat=Vi))),regimenames)

methodnames <- c("exch_lucy",
                 "indep_lucy",
                 "exch_plugin",
                 "unstr_plugin",
                 "mm_slopes",
                 "mm_intercept",
                 "trueV")

onesimrun <- function(f_sl){
  dd <- f_sl()
  d_aw <- data.frame(get_aug_weight(dd))
  d_2a <- data.frame(get_2aug(dd))
  fit_exch_lucy <- geeglm_smart_exch(d_aw, ff_fixef)
  fit_indep_lucy <- geeglm_smart_indep(d_aw, ff_fixef)
  fit_exch_plugin <- fitsmart_plugin_wr(d_aw, ff_fixef, a1s=c(1,-1), a2s=c(1,-1), corstr='exchangeable')
  fit_unstr_plugin <- fitsmart_plugin_wr(d_aw, ff_fixef, a1s=c(1,-1), a2s=c(1,-1), corstr='unstructured_a1a2')
  fit_mm_slopes <- fit_smart_lmer(d_2a, d_aw, ff_lmer_slopes, ff_fixef, ff_Z_slopes)
  fit_mm_intercept <- fit_smart_lmer(d_2a, d_aw, ff_lmer_intercept, ff_fixef, ff_Z_intercept)
  fit_trueV <- betahat_se_wr(d_aw, ff_fixef, vcomplist_a1a2_true)
  
  ## betahat
  coefmat <- 
    rbind(fit_exch_lucy$b,
          fit_indep_lucy$b,
        fit_exch_plugin$b,
        fit_unstr_plugin$b,
        fit_mm_slopes$b,
        fit_mm_intercept$b,
        fit_trueV$b)
  
  vlist_bhat <- setNames(list(fit_exch_lucy$vcov,
                         fit_indep_lucy$vcov,
                         fit_exch_plugin$vcov,
                         fit_unstr_plugin$vcov,
                         fit_mm_slopes$vcov,
                         fit_mm_intercept$vcov,
                         fit_trueV$vcov), methodnames)
  
  ## V estimates
  Vhat_mm_slopes <- get_Vhat_lmer(fit_mm_slopes, tvec, ff_Z_slopes)
  Vhat_mm_intercept <- get_Vhat_lmer(fit_mm_intercept, tvec, ff_Z_intercept)
  avg_V_unstr <- matrix(0,nrow=nrow(Vi),ncol=ncol(Vi))
  for (cregime in regimenames){
    avg_V_unstr <- avg_V_unstr + fit_unstr_plugin$Vhat_a1a2[[cregime]]
  }
  avg_V_unstr <- (1/length(regimenames)) * avg_V_unstr
  Vhat_trueV <- Vi
  ## Vhat
  vhatmat <- 
    rbind(as.numeric(fit_exch_lucy$Vhat),
          as.numeric(fit_indep_lucy$Vhat),
          as.numeric(fit_exch_plugin$Vhat_a1a2[[1]]),
          as.numeric(avg_V_unstr),
          as.numeric(Vhat_mm_slopes),
          as.numeric(Vhat_mm_intercept),
          as.numeric(Vhat_trueV))
  vnorms <-
    c(norm(fit_exch_lucy$Vhat - Vi, type='F'),
      norm(fit_indep_lucy$Vhat - Vi, type='F'),
      norm(fit_exch_plugin$Vhat_a1a2[[1]] - Vi, type='F'),
      norm(avg_V_unstr - Vi, type='F'),
      norm(Vhat_mm_slopes - Vi, type='F'),
      norm(Vhat_mm_intercept - Vi, type='F'),
      norm(Vhat_trueV - Vi, type='F'))
      
  return(list(
    bhat=coefmat,
    vlist_bhat=vlist_bhat,
    vnorms = vnorms,
    vhatmat = vhatmat,
    method=methodnames
  ))
}

simparm <- list(nsim=nsim,
                N=N,G=G, Vtruelist = Vtruelist, tvec=tvec,knot=knot,sigma=sigma,
                ff_lmer_slopes=ff_lmer_slopes, ff_fixef=ff_fixef,
                X=X, theta=theta,psi=psi,cutoff=cutoff,ff_Zgen=ff_Zgen, myseed=myseed,
                mycores=mycores,pinr_a1=pinr_a1,mmean=mmean,vardat=vardat)

f_sl_small <- datfunc_mm(N, G, tvec, knot, sigma, X, alpha_small, 
                         theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE)


f_sl_med <- datfunc_mm(N, G, tvec, knot, sigma, X, alpha_med,
                       theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE)



f_sl_large <- datfunc_mm(N, G, tvec, knot, sigma, X, alpha_large,
                         theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE)



if (doEffsize == 'small') {
    cat('\n\nsaving results for small effect size\n')
    res_small <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_small)))
    save(alpha_small, small_effsizes, truecoefs_small, res_small, simparm,
     file= paste("sim1-smalleffect-N",N,"-nsim",nsim,".RData", sep=''))
    rm(res_small)
} else if (doEffsize == 'med') {
    cat('\n\nsaving results for medium effect size\n')
    res_med <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_med)))
    save(alpha_med, med_effsizes, truecoefs_med, res_med, simparm,
     file= paste("sim1-medeffect-N",N,"-nsim",nsim,".RData", sep=''))

} else if (doEffsize == 'large') {
    cat('\n\nsaving results for large effect size\n')
    res_large <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_large)))
    save(alpha_large, large_effsizes, truecoefs_large, res_large, simparm,
     file= paste("sim1-largeeffect-N",N,"-nsim",nsim,".RData", sep=''))

} else {
    cat('\nsaving results for small, medium, and large effect sizes\n')
    res_small <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_small)))
    save(alpha_small, small_effsizes, truecoefs_small, res_small, simparm,
     file= paste("sim1-smalleffect-N",N,"-nsim",nsim,".RData", sep=''))
    rm(res_small)
    res_med <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_med)))
    save(alpha_med, med_effsizes, truecoefs_med, res_med, simparm,
     file= paste("sim1-medeffect-N",N,"-nsim",nsim,".RData", sep=''))
    rm(res_med)
    res_large <- mclapply(1:nsim, function(ix) return(onesimrun(f_sl_large)))
    save(alpha_large, large_effsizes, truecoefs_large, res_large, simparm,
     file= paste("sim1-largeeffect-N",N,"-nsim",nsim,".RData", sep=''))
}
