
get_Vi_inverse <- function(Vhatlist, cregime, tveci,uniq_times){
  V <- Vhatlist[[cregime]]
  tpos_i <- which(uniq_times %in% tveci)
  V <- V[tpos_i, tpos_i]
  return(solve(V))
}

betahat_se_wr <- function(dat_aug_weight, ff_fixef, Vhatlist, betahat = NULL, se = TRUE, Drep = NULL, Vi_inv_list=NULL) {
  if(is.null(Drep)) Drep <- model.matrix(ff_fixef, dat_aug_weight)
  p <- ncol(Drep)
  uniq_ids <- unique(dat_aug_weight$id)
  N <- length(uniq_ids)
  uniq_times <- unique(dat_aug_weight$time)
  if (!is.null(betahat)){
    muhat <- Drep %*% betahat
    dat_aug_weight$residfit <- dat_aug_weight$Y - muhat
  }
  
  if(is.null(Vi_inv_list)){
    Vi_inv_list <- setNames(vector('list', length(Vhatlist)),
                            names(Vhatlist))
  }
  
  split_dat_aug_weight <- split(dat_aug_weight,dat_aug_weight$idrep)
  cW_Vi_bdiag_list <- vector('list', length(split_dat_aug_weight))
  
  for (irep in seq_along(split_dat_aug_weight)){
    dati_regime <- split_dat_aug_weight[[irep]]
    cregime <- paste(dati_regime$A1[1],dati_regime$A2[1],sep='.')
    cW <- dati_regime$W[1]
    tveci <- dati_regime$time
    tpos_i <- which(uniq_times %in% tveci)
    tpos_i_string <- paste(tpos_i,collapse='.')
    
    if (is.null(Vi_inv_list[[cregime]][[tpos_i_string]])){
      Vi_inv_list[[cregime]][[tpos_i_string]] <- get_Vi_inverse(Vhatlist,cregime,tveci,uniq_times)
    }
    cW_Vi_bdiag_list[[irep]] <- cW * Vi_inv_list[[cregime]][[tpos_i_string]]
  }
  Drep_split_ix <- split(1:nrow(Drep),dat_aug_weight$idrep)
  
  wXi_t_Vi_inv_list <- mapply(Drep_ix = Drep_split_ix,
                              cW_Vi_inv = cW_Vi_bdiag_list,
                              function(Drep_ix, cW_Vi_inv){
                                return(crossprod(Drep[Drep_ix,], cW_Vi_inv))
                              },SIMPLIFY = FALSE)
  wXi_t_Vi_inv <- do.call('cbind', wXi_t_Vi_inv_list)
  
  A <- wXi_t_Vi_inv %*% Drep
  wXi_t_Vi_inv_Yi <- wXi_t_Vi_inv %*% dat_aug_weight$Y
  
  if (!is.null(betahat)){ # betahat provided, just return the SE
    map_idrep_id <- tapply(dat_aug_weight$idrep, dat_aug_weight$id, unique)
    Ui_a1a2_list <- 
      mapply(x1 = wXi_t_Vi_inv_list,
             Drep_ix = Drep_split_ix,
             function(x1, Drep_ix){
               return(x1 %*% dat_aug_weight$residfit[Drep_ix])
             },SIMPLIFY = F)
    
    Ui_outer_list <- 
      lapply(map_idrep_id, function(idreps){
        if(length(idreps)>1){
          return(tcrossprod(reduce(Ui_a1a2_list[idreps], `+`)))
        } else return(tcrossprod(Ui_a1a2_list[[idreps]]))
      })
    meat <- (1/N) * reduce(Ui_outer_list, `+`)
    bread <- N * solve(A)
    semat <- (1 / N) * (bread %*% meat %*% bread)
    return(semat)
  } else{ # compute betahat
    betahat <- as.numeric( solve(A, wXi_t_Vi_inv_Yi) )
    if (se){ # compute the standard errors for this betahat
      return(
        list(
          b = betahat,
          vcov = betahat_se_wr(dat_aug_weight, ff_fixef, vcomplist_a1a2, 
                               betahat=betahat, Drep = Drep, Vi_inv_list = Vi_inv_list))
      )
    } else { # just return the betahat
      return(
        list(b = betahat)
      )
    }
    
  }
}


# Compute the weighted-replicated hat(beta) coefficients
# betahat: optional pre-computed coefficients, used to compute standard errors
#          when this argument is provided, function returns standard errors only
# se: when betahat is NULL and se is FALSE, only compute hat(beta) without standard errors
#     when betahat is NULL and se = TRUE, computes hat(beta) and standard errors (requires two loops through the data)
# vcomplist_a1a2: named list of variance matrix estimates for each regime. 
#                 names correspond to regimes: "a1.a2" e.g. "1.1", "1.-1"
betahat_se_wr_old <- function(dat_aug_weight, ff_fixef, vcomplist_a1a2, betahat = NULL, se = TRUE) {
  Drep <- model.matrix(ff_fixef, dat_aug_weight)
  p <- ncol(Drep)
  uniq_ids <- unique(dat_aug_weight$id)
  N <- length(uniq_ids)
  uniq_times <- unique(dat_aug_weight$time)
  wXtVX <- matrix(0, nrow=p, ncol=p)
  wXtVY <- vector('numeric', length=p)
  if (!is.null(betahat)){
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
      solve_V <- solve(V)
      cW <- dati_regime$W[1]
      Da1a2 <- model.matrix(ff_fixef, dati_regime)
      tDa1a2 <- t(Da1a2)
      cWtDa1a2_solve_V <- cW * tDa1a2 %*% solve_V
      wXtVX <- wXtVX + cWtDa1a2_solve_V %*% Da1a2
      wXtVY <- wXtVY + cWtDa1a2_solve_V %*% Yi
      if (!is.null(betahat)){
        ri <- dati_regime$residfit
        Ui <- Ui + cW * t(Da1a2) %*% solve(V, ri)
        Ji <- Ji + cW * t(Da1a2) %*% solve_V %*% Da1a2
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
  Drep <- model.matrix(ff_fixef, dat_aug_weight)
  regimenames <- with(expand.grid(A1=a1s,A2=a2s),paste(A1,A2,sep='.'))
  maxni <- length(unique(dat_aug_weight$time))
  
  ## Start with V0=I, compute betahat0
  b0_se0 <- betahat_se_wr(dat_aug_weight, ff_fixef,
                          Vhatlist = setNames(lapply(regimenames, function(x) return(diag(maxni))),
                                              regimenames),
                         betahat = NULL, se = FALSE, Drep = Drep, Vi_inv_list=NULL)
  betahat0 <- b0_se0$b
  
  muhat0 <- as.numeric(Drep %*% betahat0) # estimated marginal mean for each replicated subject
  dat_aug_weight$resid0 <- dat_aug_weight$Y - muhat0
  p <- ncol(Drep)
  uniq_ids <- unique(dat_aug_weight$id)
  N <- length(uniq_ids)
  uniq_times <- unique(dat_aug_weight$time)
  ri0_expand <- vector('numeric', length=maxni)
  indic_tjobs <- vector('integer', length=maxni)
  nregime <- length(a1s) * length(a2s)
  vcomplist_a1a2 <- setNames(vector('list',length=nregime),
                             regimenames)
  
  Ntmat_allregime <- reduce(tapply(dat_aug_weight$time, dat_aug_weight$id,
                                function(tt) return(tcrossprod(1 * (uniq_times %in% tt)))),
                         `+`)
  Nt_allregime <- diag(Ntmat_allregime)
  ## Loop through the regimes and compute
  ## outer products of residual vectors, and sums of weights
  for (a1fix in a1s){
    for (a2fix in a2s){
      cregime <- paste(a1fix,a2fix,sep='.')
      dat_a1a2 <- subset(dat_aug_weight, A1==a1fix & A2 == a2fix)
      
      split_dat_a1a2 <- split(dat_a1a2, dat_a1a2$idrep)
      a1a2_ri_ni_W <- 
        lapply(split_dat_a1a2, function(dati_a1a2){
        ni <- nrow(dati_a1a2)
        tveci <- dati_a1a2$time
        ri0_expand[uniq_times %in% tveci] <- dati_a1a2$resid0
        ri0_expand[!(uniq_times %in% tveci)] <- 0
        indic_tjobs[uniq_times %in% tveci] <- 1
        indic_tjobs[!(uniq_times %in% tveci)] <- 0
        cW <- dati_a1a2$W[1]
        wriprodmat <- cW * tcrossprod(ri0_expand)
        indic_ts_obs_mat <- tcrossprod(indic_tjobs)
        
        return(list(
          indic_ts_obs_mat = indic_ts_obs_mat,
          nisum = ni,
          wriprodmat = wriprodmat,
          wriprodmat_ltri_scaled = wriprodmat / (sum(indic_ts_obs_mat[lower.tri(indic_ts_obs_mat,diag=F)])),
          wriprodmat_ni1_scaled = wriprodmat / (ni - 1)
        ))
      })
      
      vcomplist_a1a2[[cregime]] <- 
        reduce(a1a2_ri_ni_W, function(x1, x2){
        return(list(
          nisum = x1$ni + x2$ni,
          indic_ts_obs_mat = x1$indic_ts_obs_mat + x2$indic_ts_obs_mat,
          wriprodmat = x1$wriprodmat + x2$wriprodmat,
          wriprodmat_ltri_scaled = x1$wriprodmat_ltri_scaled + x2$wriprodmat_ltri_scaled,
          wriprodmat_ni1_scaled = x1$wriprodmat_ni1_scaled + x2$wriprodmat_ni1_scaled
        ))
      })
    }
  }
  
  vcomplist_a1a2<-
    lapply(vcomplist_a1a2, function(x){
    x$sig2_t_a1a2 <- diag(x$wriprodmat) / Nt_allregime
    x$sig2_a1a2 <- sum(diag(x$wriprodmat)) / sum(Nt_allregime)
    x$sigma_a1a2 <- sqrt(x$sig2_a1a2)
    x$sigma_t_a1a2 <- sqrt(x$sig2_t_a1a2)
    x$sigma_t_s_a1a2 <- outer(x$sigma_t_a1a2, x$sigma_t_a1a2, `*`)
    x$rho_ts_a1a2 <- x$wriprodmat[lower.tri(x$wriprodmat,diag=F)] / (Ntmat_allregime[lower.tri(Ntmat_allregime,diag=F)] *
                                                                       x$sigma_t_s_a1a2[lower.tri(x$sigma_t_s_a1a2,diag=F)])
    x$rho_a1a2 <- (1 / N) * sum(x$wriprodmat_ltri_scaled[lower.tri(x$wriprodmat,diag=F)] / (x$sigma_t_s_a1a2[lower.tri(x$sigma_t_s_a1a2,diag=F)]))
    x$psi_a1a2 <- (1 / N) * sum(x$wriprodmat_ltri_scaled[lower.tri(x$wriprodmat,diag=F)] / (x$sig2_a1a2))
    x$tau_t_a1a2 <- (1 / N) * sum(mapply(i=1:(maxni-1), j=2:maxni, function(i,j) return(x$wriprodmat_ni1_scaled[i,j] / x$sigma_t_s_a1a2[i,j])))
    x$tau_a1a2 <- (1 / N) * sum(mapply(i=1:(maxni-1), j=2:maxni, function(i,j) return(x$wriprodmat_ni1_scaled[i,j] / x$sig2_a1a2)))
    return(x)
  })
  
  ## Given the working correlation structure and those outer products,
  ## compute Vhat and the residual variance estimates
  if (corstr=='all'){
    corstr_all <- c('exchangeable_t_a1a2','exchangeable_t','exchangeable',
      'unstructured_a1a2','unstructured',
      'ar1_t_a1a2','ar1_a1a2','ar1',
      'independence_t_a1a2', 'independence_t','independence')
    Vhatlist_all <- 
      setNames(
      lapply(corstr_all,
           function(cx) return(vhat_plugin(a1s,a2s,vcomplist_a1a2, N, p, maxni,corstr=cx))
      ),
      corstr_all
    )
    ret_all <- 
      lapply(Vhatlist_all, function(Vhatlist) {
      betahat_se_final <- betahat_se_wr(dat_aug_weight, ff_fixef, Vhatlist, Drep=Drep)
      return(list(
        b=betahat_se_final$b,
        vcov=betahat_se_final$vcov,
        Vhat_a1a2=Vhatlist))
    } )
    return(ret_all)
      
  } else {
    Vhatlist <- 
      vhat_plugin(a1s, a2s, vcomplist_a1a2, N, p, maxni, corstr=corstr)
    betahat_se_final <- betahat_se_wr(dat_aug_weight, ff_fixef, Vhatlist, Drep=Drep)
    return(list(b=betahat_se_final$b,
                vcov=betahat_se_final$vcov,
                Vhat_a1a2 = Vhatlist))
  }
  
}



# Given the outer products of residuals,
# sums of weights, etc. stored in vcomplist_a1a2,
# compute the estimated Vhat for each regime given the
# desired correlation structure
vhat_plugin <- function(a1s, a2s, vcomplist_a1a2, N, p, maxni, corstr = 'exchangeable'){
  nregime <- length(vcomplist_a1a2)
  if (corstr=='exchangeable_t_a1a2'){
    ret <- 
      lapply(vcomplist_a1a2,
           function(x){
             x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
             x$Vhat[lower.tri(x$Vhat,diag=F)] <- x$sigma_t_s_a1a2[lower.tri(x$Vhat,diag=F)] * x$rho_a1a2
             diag(x$Vhat) <- x$sig2_t_a1a2
             x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
             return(x$Vhat)
           })
   return(ret) 
  } else if (corstr == 'exchangeable_t'){
    
    rho <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$rho_a1a2)),`+`)
    sig2_t <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$sig2_t_a1a2)),`+`)
    sig_t_s <- outer(sqrt(sig2_t), sqrt(sig2_t), `*`)
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- sig_t_s[lower.tri(x$Vhat, diag = F)] * rho
               diag(x$Vhat) <- sig2_t
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
  } else if (corstr=='exchangeable'){
    psi <- (1 / nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$psi_a1a2)),`+`)
    sig2 <- (1/nregime)  * reduce(lapply(vcomplist_a1a2,function(x) return(x$sig2_a1a2)),`+`)
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- sig2 * psi
               diag(x$Vhat) <- sig2
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
  } else if (corstr=='unstructured_a1a2') {
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- x$sigma_t_s_a1a2[lower.tri(x$Vhat,diag=F)] * x$rho_ts_a1a2
               diag(x$Vhat) <- x$sig2_t_a1a2
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
    return(ret) 
  }  else if (corstr=='unstructured'){
    rho_ts <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$rho_ts_a1a2)),`+`)
    sig2_t <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$sig2_t_a1a2)),`+`)
    sig_t_s <- outer(sqrt(sig2_t), sqrt(sig2_t), `*`)
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- sig_t_s[lower.tri(x$Vhat,diag=F)] * rho_ts
               diag(x$Vhat) <- sig2_t
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
    return(ret)
  } else if (corstr=='ar1_t_a1a2'){
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- x$sigma_t_s_a1a2[lower.tri(x$Vhat,diag=F)] * (x$tau_t_a1a2^abs(row(x$Vhat) - col(x$Vhat)))[lower.tri(x$Vhat,diag=F)]
               diag(x$Vhat) <- x$sig2_t_a1a2
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
    return(ret)
  } else if (corstr=='ar1_a1a2'){
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- x$sig2_a1a2 * (x$tau_a1a2^abs(row(x$Vhat) - col(x$Vhat)))[lower.tri(x$Vhat,diag=F)]
               diag(x$Vhat) <- x$sig2_a1a2
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
    return(ret)
  } else if (corstr=='ar1'){
    sig2_t <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$sig2_t_a1a2)),`+`)
    sig_t_s <- outer(sqrt(sig2_t), sqrt(sig2_t), `*`)
    tau <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$tau_a1a2)),`+`)
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               x$Vhat[lower.tri(x$Vhat,diag=F)] <- sig_t_s[lower.tri(x$Vhat,diag=F)] * (tau^abs(row(x$Vhat) - col(x$Vhat)))[lower.tri(x$Vhat,diag=F)]
               diag(x$Vhat) <- sig2_t
               x$Vhat[upper.tri(x$Vhat,diag=F)] <- x$Vhat[lower.tri(x$Vhat,diag=F)]
               return(x$Vhat)
             })
    return(ret)
  } else if (corstr=='independence_t_a1a2'){
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               diag(x$Vhat) <- x$sig2_t_a1a2
               return(x$Vhat)
             })
    return(ret)
  } else if (corstr=='independence_t') {
    sig2_t <- (1/nregime) * reduce(lapply(vcomplist_a1a2,function(x) return(x$sig2_t_a1a2)),`+`)
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               diag(x$Vhat) <- sig2_t
               return(x$Vhat)
             })
    return(ret)
  } else if (corstr=='independence') {
    sig2 <- (1/nregime)  * reduce(lapply(vcomplist_a1a2,function(x) return(x$sig2_a1a2)),`+`)
    ret <- 
      lapply(vcomplist_a1a2,
             function(x){
               x$Vhat <- matrix(0, nrow=maxni, ncol=maxni)
               diag(x$Vhat) <- sig2
               return(x$Vhat)
             })
    return(ret)
  }
  else {
    return(NULL)
  }
}
