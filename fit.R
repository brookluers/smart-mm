library(lme4)
library(geepack)

getSE_blups_lmer <- function(fitmer, ff_fixef, ff_Z, dat_aug_weight){
  betahat <- fixef(fitmer)
  p <- length(betahat)
  Drep <- model.matrix(ff_fixef, dat_aug_weight)
  muhat_i <- as.numeric(Drep %*% betahat) # estimated marginal mean for each replicated subject
  dat_aug_weight$muhat <- muhat_i
  dat_aug_weight$resid <- with(dat_aug_weight, Y - muhat)
  uniq_ids <- unique(dat_aug_weight$id)
  uniq_times <- unique(dat_aug_weight$time)
  N <- length(uniq_ids)
  J <- matrix(0, nrow=p, ncol=p)
  meat <- matrix(0, nrow=p, ncol=p)
  blups <- vector('list', length=N)
  names(blups) <- uniq_ids
  
  fitmerdims <- getME(fitmer,'devcomp')$dims
  q <- ncol(model.matrix(ff_Z, data.frame(time=uniq_times, Y=-99)))
  # Z <- getME(fitmer, "Z") 
  VC <- as.data.frame(VarCorr(fitmer), order='lower.tri')
  ltri <- subset(VC, grp!='Residual')$vcov
  Ghat <- matrix(nrow=q,ncol=q)
  Ghat[lower.tri(Ghat,diag=T)] <- ltri
  Ghat[upper.tri(Ghat)] <- Ghat[lower.tri(Ghat)]
  sig2e_hat <- subset(VC, grp=='Residual')$vcov
  
  for (i in uniq_ids){
    cdat <- subset(dat_aug_weight, id==i)  # data for subject i
    ui <- vector('numeric', p)
    Ji <- matrix(0, nrow=p, ncol=p)
    blup_numer <- 0
    blup_denom <- 0
    for (irep in unique(cdat$idrep)){ # loop over each replicated "subject" within this true subject
      cdat_regime <- subset(cdat, idrep==irep)
      ni <- nrow(cdat_regime)
      cW <- cdat_regime$W[1]
      Da1a2 <- model.matrix(ff_fixef, cdat_regime)
      # Zi <- rep(1, ni)
      # cdat_regime$Y <- -99
      Zi <- model.matrix(ff_Z, cdat_regime)
      Vi <- Zi %*% Ghat %*% t(Zi) + sig2e_hat * diag(ni)
      Viinverse <- solve(Vi)
      # blup_numer <- blup_numer + cW * t(Zi) %*% Vi %*% cdat_regime$resid
      blup_numer <- blup_numer + cW * Ghat %*% t(Zi) %*% Viinverse %*% cdat_regime$resid
      blup_denom <- blup_denom + cW
      Ji <- Ji + cW * t(Da1a2) %*% Viinverse %*% Da1a2
      ui <- ui + cW * t(Da1a2) %*% Viinverse %*% cdat_regime$resid
    }
    blups[[as.character(i)]] <- blup_numer / blup_denom
    # colnames(blups[[as.character(i)]]) <- colnames(Zi)
    J <- J + (1 / N) * Ji
    meat <- meat + (1 / N) * tcrossprod(ui)
  }
  bread <- solve(J)
  vcov <- (1 / N) *( bread %*% meat %*% bread )
  return(list(blups=blups,
              Gtri = ltri,
              sig2hat = sig2e_hat,
              vcov=vcov))
}


# ff_Z: just the random effects part of the formula, conditional on a single id
##     e.g.  Y ~ 1 for interceptso nly
# ff_lmer: the full lmer formula, including fixed and random effects
fit_smart_lmer <- function(dat_2aw, dat_aw, ff_lmer, ff_fixef, ff_Z, blups = FALSE) {
  fmer <- lmer(ff_lmer, data = dat_2aw, REML = FALSE)
  finfo <- getSE_blups_lmer(fmer, ff_fixef, ff_Z, dat_aw)
  ret <- 
    list(b=fixef(fmer),
       vcov=finfo$vcov,
       Gtri = finfo$Gtri,
       sig2hat=finfo$sig2hat)
  if (blups){
    ret$blups <- finfo$blups
  }
  return(ret)
}

get_Vhat_lmer <- function(fitmer, tvec, ff_Z){
  Gtri <- fitmer$Gtri
  ni <- length(tvec)
  Zi <- model.matrix(ff_Z, data.frame(time=tvec,Y=-99))
  q <- ncol(Zi)
  Ghat <- matrix(0, ncol=q, nrow=q)
  Ghat[lower.tri(Ghat,diag=T)] <- Gtri
  Ghat[upper.tri(Ghat,diag=F)] <- Ghat[lower.tri(Ghat,diag=F)]
  sig2hat <- fitmer$sig2hat
  return(Zi %*% Ghat %*% t(Zi) + sig2hat * diag(ni))
}

get_corrmat_exch_geem <- function(){
  require(Matrix)
  nonrepl <- matrix()
  maxni <- max(table(dat_aug_weight$id))
  cmat <- diag(maxni)
  
}

## get_zcor_exch:
##  return the custom zcor matrix to specify exchangeable correlation
## for the geeglm() function;
##    sets the correlation between artifical replicates for a single subject to zero
##    use of the default 'exchangeable' corstr in geeglm() will produce invalid esitmates + SEs
##
## arguments:
## dat_aug_weight: augmented+weighted data set
## idcol_original: string name of the column containing the original (true) subject IDs
## repl_ids: vector of unique IDs that are replicated in dat_aug_weight

# ####   documentation for the confusing genZcor function #####
# nrows = sum_i { number of upper-triangle entries in ith subject's correlation matrix } 
# ncols =   maximum number of possible upper-triangle entries for a single subject: (1/2) * ( (max_i {n_i})^2 - (max_i {n_i}) )
# ncols =   number of parameters in the unstructured correlation matrix
# each block of rows corresponds to (n_i^2 - n_i) / 2 upper-triangle entries of the ith correlation matrix
# Zcor[person i] %*% alpha gives a vector of correlation parameters for the ith subject
# myX <- genZcor(clusz=c(4,8),waves=c(1:4,1:8),corstrv=4)
# 
# Z = genZcor(...corstrv=4) for unstructured:
#      waves indicates time-ordered observations
#    
#      places 1s so that Z[1:((n_i^2 - n_i)/2), ] %*% alpha yields parameters from alpha for the ith subject
#      
get_zcor_exch <- function(dat_aug_weight, idcol_original, repl_ids){
  
  nivec_repl <- table(dat_aug_weight[[idcol_original]])
  max_nirep <- max(nivec_repl)
  uniq_ids <- unique(dat_aug_weight[[idcol_original]])
  final_zcor <- NULL
  for (cur_id in uniq_ids){
    nicur <- nivec_repl[as.character(cur_id)]
    curZcor <- 
      genZcor(clusz = nicur,
              waves = 1:nicur,
              corstrv=4)
    if (cur_id %in% repl_ids){
      half_ni <- nicur / 2
      lwrix <- half_ni
      zero_ix <- NULL
      for (ix_stepsize in (half_ni-1):0){
        uprix <- lwrix + half_ni - 1
        zero_ix <- c(zero_ix, lwrix:uprix)
        lwrix <- uprix + ix_stepsize
      }
      exch_nonzero_ix <- (1:ncol(curZcor))[-zero_ix]
      corr_design <- matrix(apply(curZcor[,exch_nonzero_ix], 1, sum),
                            nrow=nrow(curZcor),
                            byrow=T)
    } else{
      corr_design <- matrix(apply(curZcor, 1, sum),
                            nrow=nrow(curZcor),
                            byrow=T)
    }
    final_zcor <- rbind(final_zcor, corr_design)
  }
  return(final_zcor)
}

get_zcor_unstr <- function(dat_aug_weight, idcol_original, repl_ids){
  
  nivec_repl <- table(dat_aug_weight[[idcol_original]])
  max_nirep <- max(nivec_repl)
  uniq_ids <- unique(dat_aug_weight[[idcol_original]])
  zcor_list <- vector('list', length=length(uniq_ids))
  names(zcor_list) <- uniq_ids
  for (cur_id in uniq_ids){
    nicur <- nivec_repl[as.character(cur_id)]
    curZcor <- 
      genZcor(clusz = nicur,
              waves = 1:nicur,
              corstrv=4)
    if (cur_id %in% repl_ids){
      half_ni <- nicur / 2
      lwrix <- half_ni
      zero_ix <- NULL
      for (ix_stepsize in (half_ni-1):0){
        uprix <- lwrix + half_ni - 1
        zero_ix <- c(zero_ix, lwrix:uprix)
        lwrix <- uprix + ix_stepsize
      }
      nonzero_ix <- (1:ncol(curZcor))[-zero_ix]
      curZcor[,zero_ix]<-0
      corr_design <- as_tibble(curZcor)
    } else{
      corr_design <- as_tibble(curZcor)
    }
    # final_zcor <- bind_rows(final_zcor, corr_design)
    zcor_list[[as.character(cur_id)]] <- corr_design
  }
  final_zcor <- bind_rows(zcor_list) %>%
    replace(.,is.na(.), 0)
  return(as.matrix(final_zcor))
}



geeglm_smart_exch <- function(dat_aug_weight, ff_gee){
  zcor_exch <- 
    get_zcor_exch(dat_aug_weight, idcol='id', 
                  repl_ids = unique(subset(dat_aug_weight, R==1)$id))
  fgee <- geeglm(ff_gee, data = dat_aug_weight,
                 weights = W,
                 waves = waves,
                 corstr = 'userdefined',
                 id = id,
                 zcor = zcor_exch)
  sfgee <- summary(fgee)
  a1 <- setNames(sfgee$geese$corr$estimate,
                 rownames(sfgee$geese$corr))
  ni <- length(unique(dat_aug_weight$time))
  Vhat1 <- matrix(0, nrow=ni, ncol=ni)
  phihat1 <- summary(fgee)$geese$scale$estimate
  if (length(a1) > 0){
    Vhat1[lower.tri(Vhat1)] <- a1 * phihat1
    Vhat1 <- Vhat1 + t(Vhat1)
    diag(Vhat1) <- phihat1
  } else{
    diag(Vhat1) <- phihat1
  }
  
  return(
    list(
      b = sfgee$coefficients[,'Estimate'],
      vcov = sfgee$cov.unscaled,
      #se = sfgee$coefficients[,'Std.err'],
      Vhat = Vhat1,
      sig2hat = phihat1
    )
  )
}

geeglm_smart_indep <- function(dat_aug_weight, ff_gee){
  
  fgee <- 
    geeglm(ff_gee, data=dat_aug_weight,
         weights=W, waves=waves,
         corstr='independence', id=id)
  
  sfgee <- summary(fgee)
  a1 <- setNames(sfgee$geese$corr$estimate,
                 rownames(sfgee$geese$corr))
  ni <- length(unique(dat_aug_weight$time))
  Vhat1 <- matrix(0, nrow=ni, ncol=ni)
  phihat1 <- summary(fgee)$geese$scale$estimate
  if (length(a1) > 0){
    Vhat1[lower.tri(Vhat1)] <- a1 * phihat1
    Vhat1 <- Vhat1 + t(Vhat1)
    diag(Vhat1) <- phihat1
  } else{
    diag(Vhat1) <- phihat1
  }
  
  return(
    list(
      b = sfgee$coefficients[,'Estimate'],
      #se = sfgee$coefficients[,'Std.err'],
      vcov = sfgee$cov.unscaled,
      Vhat = Vhat1,
      sig2hat = phihat1
    )
  )
}

geeglm_smart_unstr <- function(dat_aug_weight, ff_gee){
  zcor_unstr <- 
    get_zcor_unstr(dat_aug_weight, idcol='id', 
                  repl_ids = unique(subset(dat_aug_weight, R==1)$id))
  fgee <- geeglm(ff_gee, data = dat_aug_weight,
                 weights = W,
                 waves = waves,
                 corstr = 'userdefined',
                 id = id,
                 zcor = zcor_unstr)
  sfgee <- summary(fgee)
  a1 <- setNames(sfgee$geese$corr$estimate,
                 rownames(sfgee$geese$corr))
  ni <- length(tvec)
  Vhat1 <- matrix(0, nrow=ni, ncol=ni)
  phihat1 <- summary(fgee)$geese$scale$estimate
  if (length(a1) > 0){
    Vhat1[lower.tri(Vhat1)] <- a1 * phihat1
    Vhat1 <- Vhat1 + t(Vhat1)
    diag(Vhat1) <- phihat1
  } else{
    diag(Vhat1) <- phihat1
  }
  
  return(
    list(
      b = sfgee$coefficients[,'Estimate'],
      #se = sfgee$coefficients[,'Std.err'],
      sfgee$cov.unscaled,
      Vhat = Vhat1,
      sig2hat = phihat1
    )
  )
}