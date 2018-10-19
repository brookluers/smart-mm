
# Compute the weighted-replicated hat(beta) coefficients
# betahat: optional pre-computed coefficients, used to compute standard errors
#          when this argument is provided, function returns standard errors only
# se: when betahat is NULL and se is FALSE, only compute hat(beta) without standard errors
#     when betahat is NULL and se = TRUE, computes hat(beta) and standard errors (requires two loops through the data)
# vcomplist_a1a2: named list of variance matrix estimates for each regime. 
#                 names correspond to regimes: "a1.a2" e.g. "1.1", "1.-1"
betahat_se_wr <- function(dat_aug_weight, ff_fixef, vcomplist_a1a2, betahat = NULL, se = TRUE) {
  Drep <- model.matrix(ff_fixef, dat_aug_weight)
  p <- ncol(Drep)
  uniq_ids <- unique(dat_aug_weight$id)
  N <- length(uniq_ids)
  uniq_times <- unique(dat_aug_weight$time)
  wXtVX <- matrix(0, nrow=p, ncol=p)
  wXtVY <- vector('numeric', length=p)
  if(!is.null(betahat)){
    muhat <- Drep %*% betahat
    dat_aug_weight$residfit <- dat_aug_weight$Y - muhat
    bread <- matrix(0, nrow=p, ncol=p)
    J <- matrix(0, nrow=p, ncol=p)
    meat <- matrix(0, nrow=p, ncol=p)
  }
  for (i in uniq_ids){
    dati <- subset(dat_aug_weight, id==i)
    if (!is.null(betahat)){
      Ji <- matrix(0, nrow=p, ncol=p)
      Ui <- vector('numeric', p)
    }
    for (irep in unique(dati$idrep)){ # loop over each replicated "subject" within this true subject
      ## This is like summing over the regimes, within a person (i) 
      dati_regime <- subset(dati, idrep==irep)
      cregime <- paste(dati_regime$A1[1], dati_regime$A2[1],sep='.')
      Yi <- dati_regime$Y
      ni <- nrow(dati_regime)
      tveci <- dati_regime$time
      V <- vcomplist_a1a2[[cregime]]$Vhat
      tpos_i <- which(uniq_times %in% tveci)
      V <- V[tpos_i, tpos_i]
      cW <- dati_regime$W[1]
      Da1a2 <- model.matrix(ff_fixef, dati_regime)
      wXtVX <- wXtVX + cW * t(Da1a2) %*% solve(V) %*% Da1a2
      wXtVY <- wXtVY + cW * t(Da1a2) %*% solve(V) %*% Yi
      if (!is.null(betahat)){
        ri <- dati_regime$residfit
        Ui <- Ui + cW * t(Da1a2) %*% solve(V, ri)
        Ji <- Ji + cW * t(Da1a2) %*% solve(V) %*% Da1a2
      }
    }
    if (!is.null(betahat)){
      J <- J + (1 / N) * Ji
      meat <- meat + (1 / N) * tcrossprod(Ui)
    }
  }
  if (!is.null(betahat)){
    bread <- solve(J)
    semat <- (1 / N) * (bread %*% meat %*% bread)
    return(semat)
  } else{
    betahat <- as.numeric(solve(wXtVX, wXtVY))
    if (se){ # compute the standard errors for this betahat
      return(
        list(
          b = betahat,
          vcov = betahat_se_wr(dat_aug_weight, ff_fixef, vcomplist_a1a2, betahat=betahat))
      )
    } else { # just return the betahat
      return(
        list(b = betahat)
      )
    }
    
  }
}


fitsmart_plugin_wr <- function(dat_aug_weight, ff_fixef, corstr='exchangeable', a1s=c(1,-1), a2s=c(1,-1)) {
  regimenames <- with(expand.grid(A1=a1s,A2=a2s),paste(A1,A2,sep='.'))
  maxni <- max(table(dat_aug_weight$idrep))
  
  ## Start with V0=I, compute betahat0
  b0_se0 <- betahat_se_wr(dat_aug_weight, ff_fixef,
                         setNames(lapply(regimenames, 
                                function(rname) {
                                  list(
                                    Vhat = diag(maxni),
                                    sig2hat = 1
                                  )
                                }), regimenames),
                         betahat = NULL, se = FALSE)
  betahat0 <- b0_se0$b
  
  Drep <- model.matrix(ff_fixef, dat_aug_weight)
  muhat0 <- as.numeric(Drep %*% betahat0) # estimated marginal mean for each replicated subject
  dat_aug_weight$resid0 <- dat_aug_weight$Y - muhat0
  p <- ncol(Drep)
  uniq_ids <- unique(dat_aug_weight$id)
  N <- length(uniq_ids)
  uniq_times <- unique(dat_aug_weight$time)
  ri0_expand <- vector('numeric', length=maxni)
  nregime <- length(a1s) * length(a2s)
  vcomplist_a1a2 <- setNames(vector('list',length=nregime),
                             regimenames)
  vcomplist_a1a2 <- 
    lapply(vcomplist_a1a2,
         function(x) return(
           list(
             sum_wi=0, # sum of Wtile_i (over i) for this regime
             sum_wni=0, # sum of Wtilde_i n_i (over i ) for this regime
             sum_wni2=0, # sum of Wtilde_i n_i^2 (over i) for this regime
             sig2hat=0, # estimated variance 
             Vhat=matrix(0, nrow=maxni, ncol=maxni),
             wriprodsum=matrix(0, nrow=maxni,ncol=maxni)
           )
         ))
  
  ## Loop through the regimes and compute
  ## outer products of residual vectors, and sums of weights
  for (a1fix in a1s){
    for (a2fix in a2s){
      cregime <- paste(a1fix,a2fix,sep='.')
      dat_a1a2 <- subset(dat_aug_weight, A1==a1fix & A2 == a2fix)
      for (irep in unique(dat_a1a2$idrep)){
        dati_a1a2 <- subset(dat_a1a2, idrep==irep)
        ni <- nrow(dati_a1a2)
        tveci <- dati_a1a2$time
        ri0_expand[uniq_times %in% tveci] <- dati_a1a2$resid0
        ri0_expand[!(uniq_times %in% tveci)] <- 0
        cW <- dati_a1a2$W[1]
        
        vcomplist_a1a2[[cregime]]$sum_wi <-
          vcomplist_a1a2[[cregime]]$sum_wi + cW
        
        vcomplist_a1a2[[cregime]]$sum_wni <- 
          vcomplist_a1a2[[cregime]]$sum_wni + cW * ni
        
        vcomplist_a1a2[[cregime]]$sum_wni2 <-
          vcomplist_a1a2[[cregime]]$sum_wni2 + cW * ni * ni
        
        vcomplist_a1a2[[cregime]]$wriprodsum <-
          vcomplist_a1a2[[cregime]]$wriprodsum + cW * tcrossprod(ri0_expand)
        
        }
    }
  }
  
  ## Given the working correlation structure and those outer products,
  ## compute Vhat and the residual variance estimates
  vcomplist_a1a2 <- 
    vhat_plugin(a1s, a2s, vcomplist_a1a2, N, p, corstr=corstr)
  betahat_se_final <- betahat_se_wr(dat_aug_weight, ff_fixef, vcomplist_a1a2)
  return(list(b=betahat_se_final$b,
              vcov=betahat_se_final$vcov,
              Vhat_a1a2 = lapply(vcomplist_a1a2, function(x) return(x$Vhat)),
              sig2hat_a1a2 = lapply(vcomplist_a1a2, function(x) return(x$sig2hat_t))))
}



# Given the outer products of residuals,
# sums of weights, etc. stored in vcomplist_a1a2,
# compute the estimated Vhat for each regime given the
# desired correlation structure
vhat_plugin <- function(a1s, a2s, vcomplist_a1a2, N, p, corstr = 'exchangeable'){
  maxni <- length(diag(vcomplist_a1a2[[1]]$Vhat))
  if (corstr=='exchangeable'){
    sig2hat_numer <- 0
    sig2hat_denom <- 0
    rhohat_numer <- 0
    rhohat_denom <- 0 
    for (a1fix in a1s){
      for (a2fix in a2s){
        cregime <- paste(a1fix, a2fix, sep='.')
        wririmat <- vcomplist_a1a2[[cregime]]$wriprodsum
        sig2hat_numer <- sig2hat_numer + sum(diag(wririmat))
        sig2hat_denom <- sig2hat_denom + vcomplist_a1a2[[cregime]]$sum_wni
        rhohat_numer <- rhohat_numer + sum(wririmat[lower.tri(wririmat,diag=F)])
        rhohat_denom <- rhohat_denom + vcomplist_a1a2[[cregime]]$sum_wni2 / 2 - 
          vcomplist_a1a2[[cregime]]$sum_wni / 2
        #vcomplist_a1a2[[cregime]]$sig2hat
        #vcomplist_a1a2[[cregime]]$Vhat
      }
    }
    sig2hat_denom <- sig2hat_denom - p
    sig2hat <- sig2hat_numer / sig2hat_denom
    rhohat_denom <- sig2hat * (rhohat_denom - p)
    rhohat <- rhohat_numer / rhohat_denom
    
    ret <- lapply(vcomplist_a1a2,
           function(x) {
             x$sig2hat <- sig2hat
             x$sig2hat_t <- rep(sig2hat, maxni)
             x$Vhat[lower.tri(x$Vhat)] <- rhohat
             x$Vhat[upper.tri(x$Vhat)] <- rhohat
             diag(x$Vhat) <- 1
             x$Vhat <- x$Vhat * sig2hat
             return(x)
           }
             )
   return(ret) 
  } else if (corstr=='unstructured_a1a2') {
    # for each regime,
    #   estimate the variance at each time point
    #   estimate the correlation individually for each pair of time points (t, s)
    sig2hat_t_numer <- vector('numeric', maxni)
    sig2hat_denom <- 0
    rhohat_ts_numer <- vector('numeric', maxni * (maxni - 1) / 2)
    rhohat_ts_denom <- vector('numeric', maxni * (maxni - 1) /2 )
    
    ret <-
      lapply(vcomplist_a1a2, function(x){
        wririmat <- x$wriprodsum
        x$sig2hat_t <- diag(wririmat) / (x$sum_wi - p)
        sighat_t <- sqrt(x$sig2hat_t)
        rho_ts_numer <- as.numeric(wririmat[lower.tri(wririmat, diag=F)])
        sig_t_sig_s <- outer(sighat_t, sighat_t,`*`)[lower.tri(x$Vhat,diag=F)]
        rho_ts_denom <- N * sig_t_sig_s
        x$rhohat_ts <- rho_ts_numer / rho_ts_denom
        x$Vhat[lower.tri(x$Vhat,diag=F)] <- x$rhohat_ts * sig_t_sig_s
        x$Vhat <- x$Vhat + t(x$Vhat)
        diag(x$Vhat) <- x$sig2hat_t
        return(x)
      })
    
    return(ret) 
  } else {
    return(NULL)
  }
}
