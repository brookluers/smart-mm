source("gen.R")
source("fit.R")
source("fitplugin.R")
library(parallel)

get_simparm <- function(cmd_args, a1s, a2s, alphalist, effsizenames,
                        psi, theta, tvec, knot, sigma, cutoff, 
                        ff_Zgen, G, covfunc_epsilon, corstr='all', missFunc = NULL){
  N <- as.numeric(cmd_args[1])
  nsim <- as.numeric(cmd_args[2])
  myseed <- as.numeric(cmd_args[3])
  mycores <- cmd_args[4]
  doEffsize <- cmd_args[5]
  set.seed(myseed)
  cat("\nN = "); cat(N)
  cat("\nnsim = "); cat(nsim)
  if (is.na(doEffsize)) {
    doEffsize <- 'all'
  }
  cat("\neffect size = "); cat(doEffsize)
  if (is.null(missFunc)){
    cat("\nNo missingness\n")
  } else {
    cat("\nmissingness function supplied\n")
  }
  
  cat("\nseed = "); cat(myseed)
  cat("\ndesired cores = "); cat(mycores);
  cat("\ndetectCores() = "); cat(detectCores())
  mycores <- max(c(mycores, detectCores() - 1))
  options(cores=mycores)
  options(mc.cores=mycores)
  cat("\n--using "); cat(mycores); cat(" cores")
  
  allregime <- expand.grid(a1 = a1s, a2 = a2s)
  regimenames <- as.character(with(allregime, interaction(a1,a2)))
  
  X <- vector('numeric', length=N)
  X[1:floor(N/2)] <- 1
  X[(floor(N/2)+1):N] <- -1
  Zi <- model.matrix(ff_Zgen, data.frame(time=tvec,Y=-99))
  ni <- length(tvec)
  regimes_effnames <- expand.grid(effsize = effsizenames,
                                  a1=a1s, a2=a2s)
  Vtruelist <- 
    lapply(effsizenames, function(effsize) {
      ret <-
        mapply(a1=allregime$a1, a2= allregime$a2,
               FUN = function(a1, a2){
                 SigmaTrueX1 <- matrix(0, nrow=ni, ncol=ni)
                 SigmaTrueX1 <- 
                   matrix(mapply(t_ix = row(SigmaTrueX1),
                                 s_ix = col(SigmaTrueX1),
                                 FUN = function(t_ix, s_ix) {
                                   cij <- get_margcov(tval = tvec[t_ix], sval=tvec[s_ix],
                                                      alpha = alphalist[[effsize]], 
                                                      psi=psi, knot=knot, a1=a1, a2=a2,
                                                      G=G, ff_Z=ff_Zgen, sigma=sigma,cutoff=cutoff,
                                                      fn_tscov = covfunc_epsilon)
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
  
  mmeanlist <- 
    lapply(effsizenames, function(effsize){
      return(
        get_mmean(alphalist[[effsize]], 
                unique(X), theta, knot, tvec, a1 = a1s, a2 = a2s, 
                G, ff_Zgen, sigma, cutoff)
      )
    })
    
  truecoeflist <- lapply(alphalist, function(aa){
    return(get_truecoefs(aa, theta, knot, G, ff_Zgen, sigma, cutoff))
  })
  
   pinr_a1 <-
      lapply(effsizenames, function(effsize){
        setNames(
          sapply(a1s, function(aa11) return(
          get_pinr(alphalist[[effsize]], knot, a1 = aa11, G, ff_Zgen, sigma, cutoff)
        )), as.character(a1s))
      })
  
  names(Vtruelist) <- 
    names(mmeanlist) <- 
    names(pinr_a1) <- effsizenames
  
  
  vardat <- expand.grid(tval=tvec,
                        a1 = a1s,
                        a2 = a2s,
                        effsize = effsizenames)
  vardat$mvar <- 
    mapply(tval = vardat$tval,
           a1=vardat$a1,
           a2=vardat$a2,
           effsize = vardat$effsize,
           FUN = function(tval, a1, a2, effsize) return(get_margvar(tval =tval, 
                                                                    alphalist[[effsize]], 
                                                                    psi, knot, 
                                                                    a1=a1, a2=a2, G, ff_Zgen, sigma,cutoff,fn_tscov = covfunc_epsilon)))
  
  Sigma_epsilon <- get_Sigma_eps(tvec, fn_tscov = covfunc_epsilon, sigma=sigma)
  
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
  
  effsizelist <-
    setNames(
      lapply(effsizenames, function(effsize){ # for each effect size
      return(
        apply(combn(1:nrow(allregime),2), 2, function(ij) { # for each regime pair
          regimemat <- as.matrix(allregime)
          a11 <- regimemat[ij[1], 1]
          a21 <- regimemat[ij[1], 2]
          a12 <- regimemat[ij[2], 1]
          a22 <- regimemat[ij[2], 2]
          return(
            get_effsizes(a11=a11, a21=a21, a12=a12, a22=a22, alpha= alphalist[[effsize]], 
                       psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                       fn_tscov = covfunc_epsilon))
        })
      )
    }), effsizenames)
  
  methodnames <- c("mm_slopes", "mm_intercept", "trueV")
  
  fgenlist <- 
    setNames(
      lapply(effsizenames, function(effsize){
    return(
      datfunc_mm(N, G, tvec, knot, sigma, X, alphalist[[effsize]], 
               theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE,
               Sigma_epsilon)
    )
  }), effsizenames)
  
  simparm <- list(
    nsim = nsim,
    effsizelist = effsizelist,
    doEffsize = doEffsize,
    corstr = corstr,
    fgenlist = fgenlist,
    missFunc = missFunc,
    N = N,
    G = G,
    Vtruelist = Vtruelist,
    tvec = tvec,
    knot = knot,
    sigma = sigma,
    ff_lmer_slopes = ff_lmer_slopes,
    ff_fixef = ff_fixef,
    alphalist = alphalist,
    X = X,
    theta = theta,
    truecoefs = truecoeflist,
    psi = psi,
    cutoff = cutoff,
    ff_Zgen = ff_Zgen, ff_Z_intercept = ff_Z_intercept,
    ff_Z_slopes = ff_Z_slopes,
    ff_lmer_intercept = ff_lmer_intercept,
    ff_lmer_slopes = ff_lmer_slopes,
    myseed = myseed,
    mycores = mycores,
    methodnames = methodnames,
    regimenames = regimenames,
    effsizenames = effsizenames,
    pinr_a1 = pinr_a1,
    mmeanlist = mmeanlist,
    a1s = a1s, a2s = a2s,
    vardat = vardat,
    covfunc_epsilon = covfunc_epsilon,
    Sigma_epsilon = Sigma_epsilon
  )
  
  return(simparm)
}

get_onesimrun <- function(simparm, effsize){
  if (!(effsize %in% simparm$effsizenames)){
    stop("incorrect effect size name passed to get_onesimrun\n")
  }
  f_sl <- simparm$fgenlist[[effsize]]
  Vtrue_a1a2 <- simparm$Vtruelist[[effsize]]
  
  onesimrun <- function(){
    dd <- f_sl()
    if (!is.null(simparm$missFunc)){
      dd <- simparm$missFunc(dd)
    }
    d_aw <- get_aug_weight(dd)
    d_2a <- data.frame(get_2aug(d_aw))
    d_aw <- data.frame(d_aw)
    
    fit_all_plugin <- fitsmart_plugin_wr(d_aw, simparm$ff_fixef, 
                                         corstr = simparm$corstr, 
                                         a1s = simparm$a1s, 
                                         a2s = simparm$a2s)
    
    fit_mm_slopes <- fit_smart_lmer(d_2a, d_aw, simparm$ff_lmer_slopes, simparm$ff_fixef, simparm$ff_Z_slopes)
    fit_mm_intercept <- fit_smart_lmer(d_2a, d_aw, simparm$ff_lmer_intercept, simparm$ff_fixef, simparm$ff_Z_intercept)
    fit_trueV <- betahat_se_wr(d_aw, simparm$ff_fixef, Vtrue_a1a2)
    
    ## betahat
    coefmat <- 
      rbind('mm_slopes' = fit_mm_slopes$b,
            'mm_intercept' = fit_mm_intercept$b,
            'trueV' = fit_trueV$b,
            do.call('rbind', lapply(fit_all_plugin, function(fit_j) return(fit_j$b))))
    # rownames(coefmat) <- simparm$methodnames
    
    vlist_bhat <- c(list('mm_slopes' = fit_mm_slopes$vcov,
                          'mm_intercept' = fit_mm_intercept$vcov,
                            'trueV' = fit_trueV$vcov),
                             lapply(fit_all_plugin, function(fit_j) return(fit_j$vcov)))
    
    ## V estimates
    Vhat_mm_slopes <- get_Vhat_lmer(fit_mm_slopes, simparm$tvec, simparm$ff_Z_slopes)
    Vhat_mm_intercept <- get_Vhat_lmer(fit_mm_intercept, simparm$tvec, simparm$ff_Z_intercept)
    
    
    Vhat_a1a2_allfits <-
      c(list('mm_slopes' = setNames(lapply(simparm$regimenames, function(cr) return(Vhat_mm_slopes)), simparm$regimenames),
                      'mm_intercept' = setNames(lapply(simparm$regimenames, function(cr) return(Vhat_mm_intercept)), simparm$regimenames),
                      'trueV' = Vtrue_a1a2),
                 lapply(fit_all_plugin, function(fj) return(fj$Vhat_a1a2)))
    avg_Vnorms <- 
      sapply(Vhat_a1a2_allfits,
             function(vhl_a1a2) {
               vnorms_a1a2 <- 
                 mapply(vhat = vhl_a1a2,
                        vtrue = Vtrue_a1a2,
                        function(vhat, vtrue) return(norm(vhat - vtrue,type='F')),
                        SIMPLIFY = F)
               return((1/length(simparm$regimenames)) * reduce(vnorms_a1a2, `+`))
             })
    
    
    vhatmat <- setNames(lapply(simparm$regimenames, function(rr)
      return(do.call(
        'rbind', lapply(Vhat_a1a2_allfits, function(vhj)
          return(vhj[[rr]][lower.tri(vhj[[rr]], diag = T)]))
      ))), simparm$regimenames)
    
    return(list(
      bhat = coefmat,
      vlist = vlist_bhat,
      vhatmat = vhatmat,
      vnorms = avg_Vnorms,
      method = rownames(coefmat)
    ))
  }
  return(onesimrun)

}


runsim <- function(simparm, simname){
  doEffsize <- simparm$doEffsize
  
  onesimrun_list <- 
    setNames(
      lapply(simparm$effsizenames, function(effsize){
    return(get_onesimrun(simparm, effsize = effsize))
  }), simparm$effsizenames)
  
  if (doEffsize == 'all'){
    doEffsize <- simparm$effsizenames
  }
  res <- setNames(vector('list', length(doEffsize)),
                  doEffsize)
  
  for (effname in doEffsize){
    cat(paste('\n\nsaving results for ', effname, ' effect size\n',sep=''))
    res <- mclapply(1:simparm$nsim, function(i) return(onesimrun_list[[effname]]()))
    save(res, simparm, 
         file= paste(simname, "-", effname, "effect-N", 
                     simparm$N, "-nsim", simparm$nsim, ".RData", sep=''))
  }
}