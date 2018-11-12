allregime <- expand.grid(a1=a1s,a2=a2s)
regimenames <- as.character(with(allregime,interaction(a1,a2)))
coefnames <- c(paste('beta', 0:6,sep=''), 'eta')

X <- vector('numeric', length=N)
X[1:floor(N/2)] <- 1
X[(floor(N/2)+1):N] <- -1

Zi <- model.matrix(ff_Zgen, data.frame(time=tvec,Y=-99))

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
names(Vtruelist) <- c('small','med','large')

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
                      a2=c(1,-1),
                      effsize = c('small','med','large'))
vardat$mvar <- 
  mapply(tval = vardat$tval,
         a1=vardat$a1,
         a2=vardat$a2,
         effsize = vardat$effsize,
         FUN = function(tval, a1, a2, effsize) return(get_margvar(tval =tval, 
                                                         get(paste('alpha_',effsize,sep='')), 
                                                         psi, knot, 
                                                         a1=a1, a2=a2, G, ff_Zgen, sigma,cutoff,fn_tscov = covfunc_epsilon)))

cat("Covariance of residual errors:\n")
( Sigma_epsilon <- get_Sigma_eps(tvec, fn_tscov = covfunc_epsilon, sigma=sigma) )

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



small_effsizes <- 
  rbind(
    get_effsizes(a11=1, a21=1, a12=-1, a22=1, alpha= alpha_small, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=1, a21=1, a12=-1, a22=-1, alpha= alpha_small, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=1, a21=-1, a12=-1, a22=1, alpha= alpha_small, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=-1, a21=-1, a12=-1, a22=1, alpha= alpha_small, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon))


med_effsizes <- 
  rbind(
    get_effsizes(a11=1, a21=1, a12=-1, a22=1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=1, a21=1, a12=-1, a22=-1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=1, a21=-1, a12=-1, a22=1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=-1, a21=-1, a12=-1, a22=1, alpha= alpha_med, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon))



large_effsizes <- 
  rbind(
    get_effsizes(a11=1, a21=1, a12=-1, a22=1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=1, a21=1, a12=-1, a22=-1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=1, a21=-1, a12=-1, a22=1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon),
    get_effsizes(a11=-1, a21=-1, a12=-1, a22=1, alpha= alpha_large, 
                 psi=psi, X=1, theta, knot, tvec=max(tvec), G, ff_Zgen, sigma=sigma, cutoff,
                 fn_tscov = covfunc_epsilon))


methodnames <- c(# "exch_lucy",
                 "indep_lucy",
                 "exch_plugin",
                 "unstr_plugin",
                 "mm_slopes",
                 "mm_intercept",
                 "trueV")

onesimrun <- function(f_sl, Vtrue_a1a2, missFunc = NULL){
  vcomplist_a1a2_true <- lapply(Vtrue_a1a2, function(x) return(list(Vhat=x)))
  
  dd <- f_sl()
  if (!is.null(missFunc)){
    dd <- missFunc(dd)
  }
  d_aw <- get_aug_weight(dd)
  d_2a <- data.frame(get_2aug(d_aw))
  d_aw <- data.frame(d_aw)
  # fit_exch_lucy <- geeglm_smart_exch(d_aw, ff_fixef)
  fit_indep_lucy <- geeglm_smart_indep(d_aw, ff_fixef)
  fit_exch_plugin <- fitsmart_plugin_wr(d_aw, ff_fixef, a1s=c(1,-1), a2s=c(1,-1), corstr='exchangeable')
  fit_unstr_plugin <- fitsmart_plugin_wr(d_aw, ff_fixef, a1s=c(1,-1), a2s=c(1,-1), corstr='unstructured_a1a2')
  fit_mm_slopes <- fit_smart_lmer(d_2a, d_aw, ff_lmer_slopes, ff_fixef, ff_Z_slopes)
  fit_mm_intercept <- fit_smart_lmer(d_2a, d_aw, ff_lmer_intercept, ff_fixef, ff_Z_intercept)
  fit_trueV <- betahat_se_wr(d_aw, ff_fixef, vcomplist_a1a2_true)
  
  ## betahat
  coefmat <- 
    rbind(# fit_exch_lucy$b,
          fit_indep_lucy$b,
          fit_exch_plugin$b,
          fit_unstr_plugin$b,
          fit_mm_slopes$b,
          fit_mm_intercept$b,
          fit_trueV$b)
  
  vlist_bhat <- setNames(list(# fit_exch_lucy$vcov,
                              fit_indep_lucy$vcov,
                              fit_exch_plugin$vcov,
                              fit_unstr_plugin$vcov,
                              fit_mm_slopes$vcov,
                              fit_mm_intercept$vcov,
                              fit_trueV$vcov), methodnames)
  
  ## V estimates
  Vhat_mm_slopes <- get_Vhat_lmer(fit_mm_slopes, tvec, ff_Z_slopes)
  Vhat_mm_intercept <- get_Vhat_lmer(fit_mm_intercept, tvec, ff_Z_intercept)
  avg_Vtrue <- matrix(0,nrow=nrow(Vhat_mm_slopes), ncol=ncol(Vhat_mm_slopes))
  avg_V_unstr <- matrix(0,nrow=nrow(Vhat_mm_slopes), ncol=ncol(Vhat_mm_slopes))
  for (cregime in regimenames){
    avg_V_unstr <- avg_V_unstr + fit_unstr_plugin$Vhat_a1a2[[cregime]]
    avg_Vtrue <- avg_Vtrue + Vtrue_a1a2[[cregime]]
  }
  avg_Vtrue <- (1 / length(regimenames)) * avg_Vtrue
  avg_V_unstr <- (1/length(regimenames)) * avg_V_unstr
  avg_Vnorms <- setNames(vector('numeric', length(methodnames)), methodnames)
  for (cregime in regimenames){
    avg_Vnorms['unstr_plugin'] <- avg_Vnorms['unstr_plugin'] + 
      (1 / length(regimenames)) * norm(fit_unstr_plugin$Vhat_a1a2[[cregime]] - Vtrue_a1a2[[cregime]],type='F')
   # avg_Vnorms['exch_lucy'] <- avg_Vnorms['exch_lucy'] + 
      #(1 / length(regimenames)) * norm(fit_exch_lucy$Vhat - Vtrue_a1a2[[cregime]], type='F')
    avg_Vnorms['indep_lucy'] <- avg_Vnorms['indep_lucy'] + 
      (1 / length(regimenames)) * norm(fit_indep_lucy$Vhat - Vtrue_a1a2[[cregime]], type='F')
    avg_Vnorms['mm_slopes'] <- avg_Vnorms['mm_slopes'] + 
      (1 / length(regimenames)) * norm(Vhat_mm_slopes - Vtrue_a1a2[[cregime]], type='F')
    avg_Vnorms['mm_intercept'] <- avg_Vnorms['mm_intercept'] + 
      (1 / length(regimenames)) * norm(Vhat_mm_intercept - Vtrue_a1a2[[cregime]], type='F')
    avg_Vnorms['exch_plugin'] <- avg_Vnorms['exch_plugin'] +
      (1 / length(regimenames)) * norm(fit_exch_plugin$Vhat_a1a2[[cregime]] - Vtrue_a1a2[[cregime]], type='F')
    avg_Vnorms['trueV'] <- 0
  }
  ## Vhat
  vhatmat <- 
    rbind(#as.numeric(fit_exch_lucy$Vhat),
          as.numeric(fit_indep_lucy$Vhat),
          as.numeric(fit_exch_plugin$Vhat_a1a2[[1]]),
          as.numeric(avg_V_unstr),
          as.numeric(Vhat_mm_slopes),
          as.numeric(Vhat_mm_intercept),
          as.numeric(avg_Vtrue))
  
  return(list(
    bhat=coefmat,
    vlist=vlist_bhat,
    vhatmat = vhatmat,
    vnorms = avg_Vnorms,
    method=methodnames
  ))
}

simparm <- list(nsim=nsim,
                N=N,G=G, Vtruelist=Vtruelist, tvec=tvec,knot=knot,sigma=sigma,ff_lmer_slopes=ff_lmer_slopes,ff_fixef=ff_fixef,
                X=X, theta=theta,psi=psi,cutoff=cutoff,ff_Zgen=ff_Zgen,myseed=myseed,
                mycores=mycores,pinr_a1=pinr_a1,mmean=mmean,vardat=vardat,
                covfunc_epsilon=covfunc_epsilon, Sigma_epsilon=Sigma_epsilon)
if (exists('missprob')){
  simparm <- c(simparm, missprob=missprob)
}
if (exists('missFunc')){
  simparm <- c(simparm, missFunc=missFunc)
}
f_sl_small <- datfunc_mm(N, G, tvec, knot, sigma, X, alpha_small, 
                         theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE,
                         Sigma_epsilon)

f_sl_med <- datfunc_mm(N, G, tvec, knot, sigma, X, alpha_med,
                       theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE,
                       Sigma_epsilon)

f_sl_large <- datfunc_mm(N, G, tvec, knot, sigma, X, alpha_large,
                         theta, psi, cutoff, ff_Z=ff_Zgen, return_po = FALSE,
                         Sigma_epsilon)

