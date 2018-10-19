library(MASS)

get_pinr_old <- function(alpha, knot, a1, G, sigma, cutoff){
  mu <- alpha[1] + knot * alpha[2] + knot * alpha[3] * a1
  zt <- c(1, knot)
  v <- as.numeric(t(zt) %*% G %*% zt) + sigma^2
  return(pnorm((cutoff - mu) / sqrt(v)))
}

get_pinr <- function(alpha, knot, a1, G, ff_Z, sigma, cutoff) {
  mu <- alpha[1] + knot * alpha[2] + knot * alpha[3] * a1
  zknotT <- model.matrix(ff_Z, data.frame(time=knot,Y=-99))
  v <- as.numeric(zknotT %*% G %*% t(zknotT) + sigma^2)
  return(pnorm((cutoff - mu) / sqrt(v)))
}

get_cutoff <- function(alpha, knot, a1, G, ff_Z, sigma, pinr) {
  mu <- alpha[1] + knot * alpha[2] + knot * alpha[3] * a1
  zknotT <- model.matrix(ff_Z, data.frame(time=knot,Y=-99))
  v <- as.numeric(zknotT %*% G %*% t(zknotT) + sigma^2)
  return(qnorm(pinr, mean=mu, sd=sqrt(v)))
}
# 
# get_cutoff <- function(alpha, knot, a1, G, sigma, pinr){
#   mu <- alpha[1] + knot * alpha[2] + knot * alpha[3] * a1
#   zt <- c(1, knot)
#   v <- t(zt) %*% G %*% zt
#   return(qnorm(pinr, mean=mu, sd=sqrt(v)))
# }

get_margvar <- function(tval, alpha, psi, knot, a1, a2, G, ff_Z, sigma, cutoff) {
  return(get_margcov(tval, tval, alpha, psi, knot, a1, a2, G, ff_Z, sigma, cutoff))
}

###
# get_margvar <- function(tval, alpha, psi, knot, a1, a2, G, sigma, cutoff){
#   require(tmvtnorm)
#   pinr <- get_pinr(alpha, knot, a1, G, sigma, cutoff)
#   if (tval <= knot){
#     return(G[1,1] + 2 * tval * G[1,2] + tval^2 * G[2,2] + sigma^2)
#   } else{
#     WWvar <- 
#       matrix(c(
#         G[1,1] + 2 * tval * G[1,2] + tval^2 * G[2,2] + sigma^2,
#         G[1,1] + knot * G[1,2] + tval * G[1,2] + tval * knot * G[2,2],
#         G[1,1] + knot * G[1,2] + tval * G[1,2] + tval * knot * G[2,2],
#         G[1,1] + 2 * knot * G[1,2] + knot^2 * G[2,2] + sigma^2),
#         nrow=2,byrow=T)
#     C1 <- (tval - knot) * alpha[6] * a2
#     C2 <- (tval - knot) * alpha[7] * a1 * a2
#     C3 <- (tval - knot) * 1 * (a1==1) * psi['1'] + 1 * (a1==-1) * psi['-1']
#     
#     eW1W2big <- 
#       mtmvnorm(mean = c(0, 0),
#                sigma = WWvar,
#                lower = c(-Inf, cutoff - alpha[1] - knot * (alpha[2] + alpha[3] * a1)),
#                upper = c(Inf, Inf))$tmean[1]
#     return(
#       G[1,1] + tval^2 * G[2,2] + 2 * tval * G[1,2] + sigma^2 +
#         (C1 + C2 - C3)^2 * pinr * (1 - pinr) -
#         2 * (C1 + C2 - C3) * (1 - pinr) * eW1W2big
#     )
#       
#   }
# }
###

get_margcov <- function(tval, sval, alpha, psi, knot, a1, a2, G, ff_Z, sigma, cutoff) {
  pinr <- get_pinr(alpha, knot, a1, G, ff_Z, sigma, cutoff)
  ztT <- model.matrix(ff_Z, data.frame(Y=-99,time=tval))
  zsT <- model.matrix(ff_Z, data.frame(Y=-99,time=sval))
  
  Ct1 <- 1 * (tval > knot) * (tval - knot) * alpha[6] * a2
  Ct2 <- 1 * (tval > knot ) * (tval - knot) * alpha[7] * a1 * a2
  Ct3 <- 1 * (tval > knot ) * (tval - knot) * 1 * (a1==1) * psi['1'] + 1 * (a1==-1) * psi['-1']
  Cs1 <- 1 * (sval > knot) * (sval - knot) * alpha[6] * a2
  Cs2 <- 1 * (sval > knot ) * (sval - knot) * alpha[7] * a1 * a2
  Cs3 <- 1 * (sval > knot ) * (sval - knot) * 1 * (a1==1) * psi['1'] + 1 * (a1==-1) * psi['-1']
  
  Wtmean_R1 <- tmoments_WtWknot(tval,alpha,knot,a1,a2,G,ff_Z,sigma,cutoff,Req1=TRUE)$tmean[1]
  Wsmean_R1 <- tmoments_WtWknot(sval,alpha,knot,a1,a2,G,ff_Z,sigma,cutoff,Req1=TRUE)$tmean[1]
  
  cov_WtR <-  (1 - pinr) * Wtmean_R1
  cov_WsR <-  (1 - pinr) * Wsmean_R1
  
  return( ztT %*% G %*% t(zsT) + 1* (sval==tval) * sigma^2 -
            (Cs1 + Cs2 - Cs3) * cov_WtR -
            (Ct1 + Ct2 - Ct3) * cov_WsR +
            (Cs1 + Cs2 - Cs3) * (Ct1 + Ct2 - Ct3) * pinr * (1 - pinr)
  )
}

get_margcov_old <- function(tval, sval, alpha, theta, psi, knot, Xval, a1, a2, G, sigma, cutoff) {
  pinr <- get_pinr_old(alpha, knot, a1, G, sigma, cutoff)
  Ct1 <- 1 * (tval > knot) * (tval - knot) * alpha[6] * a2
  Ct2 <- 1 * (tval > knot ) * (tval - knot) * alpha[7] * a1 * a2
  Ct3 <- 1 * (tval > knot ) * (tval - knot) * 1 * (a1==1) * psi['1'] + 1 * (a1==-1) * psi['-1']
  Cs1 <- 1 * (sval > knot) * (sval - knot) * alpha[6] * a2
  Cs2 <- 1 * (sval > knot ) * (sval - knot) * alpha[7] * a1 * a2
  Cs3 <- 1 * (sval > knot ) * (sval - knot) * 1 * (a1==1) * psi['1'] + 1 * (a1==-1) * psi['-1']

  Wtmean_R1 <- tmoments_WtWknot_old(tval,alpha,knot,a1,a2,G,sigma,cutoff,Req1=TRUE)$tmean[1]
  Wsmean_R1 <- tmoments_WtWknot_old(sval,alpha,knot,a1,a2,G,sigma,cutoff,Req1=TRUE)$tmean[1]

  cov_WtR <-  (1 - pinr) * Wtmean_R1
  cov_WsR <-  (1 - pinr) * Wsmean_R1

  return(
    G[1,1] + (tval + sval) * G[1,2] + tval*sval * G[2,2]  + 1 * (sval == tval) * sigma^2 -
      (Cs1 + Cs2 - Cs3) * cov_WtR -
      (Ct1 + Ct2 - Ct3) * cov_WsR +
      (Cs1 + Cs2 - Cs3) * (Ct1 + Ct2 - Ct3) * pinr * (1 - pinr)
  )
}

get_vcovmat <- function(a1, a2, tvec, alpha, theta, psi, knot, G, ff_Z, sigma, cutoff) {
  ni <- length(tvec)
  a1 <-1; a2 <- 1
  SigmaTrueX1 <- matrix(0, nrow=ni, ncol=ni)
  SigmaTrueX1 <- 
    matrix(mapply(t_ix = row(SigmaTrueX1),
                  s_ix = col(SigmaTrueX1),
                  FUN = function(t_ix, s_ix) {
                    cij <- get_margcov(tval=tvec[t_ix], sval=tvec[s_ix],
                                       alpha,psi,knot,a1,a2,
                                       G,ff_Z, sigma,cutoff)
                    return(
                      cij
                    )
                  }),
           nrow=ni,byrow=T)
  return(SigmaTrueX1)
}

tmoments_WtWknot <- function(tval, alpha, knot, a1, a2, G, ff_Z, sigma, cutoff, Req1=TRUE){
  require(tmvtnorm)
  ZtT <- model.matrix(ff_Z, data.frame(Y = -99, time=tval))
  ZknotT <- model.matrix(ff_Z, data.frame(Y=-99, time=knot))
  
  WWvar <- 
    matrix(
      c(ZtT %*% G %*% t(ZtT) + sigma^2, ZknotT %*% G %*% t(ZtT),
        ZknotT %*% G %*% t(ZtT), ZknotT %*% G %*% t(ZknotT) + sigma^2),
      nrow=2,byrow=T)
  
  if (Req1){
    lower <- c(-Inf, cutoff - alpha[1] - knot * (alpha[2] + alpha[3] * a1))
    upper <- c(Inf, Inf)
  } else{
    lower <- c(-Inf, -Inf)
    upper <- c(Inf, cutoff - alpha[1] - knot*(alpha[2] + alpha[3] * a1))
  }
  tmoments <- 
    mtmvnorm(mean = c(0, 0),
             sigma = WWvar,
             lower = lower,
             upper = upper)
  return(tmoments)
  
}

tmoments_WtWknot_old <- function(tval, alpha, knot, a1, a2, G, sigma, cutoff, Req1 = TRUE) {
  require(tmvtnorm)
  WWvar <-
    matrix(c(
      G[1,1] + 2 * tval * G[1,2] + tval^2 * G[2,2] + sigma^2,
      G[1,1] + knot * G[1,2] + tval * G[1,2] + tval * knot * G[2,2],
      G[1,1] + knot * G[1,2] + tval * G[1,2] + tval * knot * G[2,2],
      G[1,1] + 2 * knot * G[1,2] + knot^2 * G[2,2] + sigma^2),
      nrow=2,byrow=T)
  if (Req1){
    lower <- c(-Inf, cutoff - alpha[1] - knot * (alpha[2] + alpha[3] * a1))
    upper <- c(Inf, Inf)
  } else{
    lower <- c(-Inf, -Inf)
    upper <- c(Inf, cutoff - alpha[1] - knot*(alpha[2] + alpha[3] * a1))
  }
  tmoments <-
    mtmvnorm(mean = c(0, 0),
             sigma = WWvar,
             lower = lower,
             upper = upper)
  return(tmoments)
}

get_mmean <- function(alpha, X, theta, knot, tvec, a1, a2, G, ff_Z, sigma, cutoff) {
  res <- as.matrix(expand.grid(a1=a1,a2=a2,time=tvec, X=X))
  mu <- mapply(a1=res[,'a1'], a2=res[,'a2'], time=res[,'time'],X=res[,'X'],
         function(a1,a2,time,X){
           pinr <- get_pinr(alpha, knot, a1, G, ff_Z, sigma, cutoff)
           alpha[1] + 1 * (time <= knot) * time * (alpha[2] + alpha[3] * a1) +
             1 * (time > knot) * knot * (alpha[2] + alpha[3] * a1) +
             1 * (time > knot) * (time - knot) * (alpha[4] + alpha[5] * a1 + alpha[6] * pinr * a2 +
                                                    alpha[7] * pinr * a1 * a2) +
             theta * X
         },
         SIMPLIFY = TRUE)
  return(cbind(res,mu=mu))
}

get_truecoefs <- function(alpha, theta, knot, G, ff_Z, sigma, cutoff) {
  pinr_1 <- get_pinr(alpha, knot, a1=1, G, ff_Z, sigma, cutoff)
  pinr_m1 <- get_pinr(alpha, knot, a1=-1, G, ff_Z, sigma, cutoff)
  
  ret <- vector('numeric', length=length(alpha)+1)
  ret[1:5] <- alpha[1:5]
  ret[6] <- alpha[6] * (pinr_1 / 2 + pinr_m1 / 2) + alpha[7] * (pinr_1/2 - pinr_m1/2)
  ret[7] <- alpha[6] * (pinr_1 / 2 - pinr_m1/2) + alpha[7] * (pinr_1/2 + pinr_m1/2)
  ret[8] <- theta
  names(ret) <- c(paste('beta',0:6,sep=''), 'eta')
  return(ret)
}

get_effsizes <- function(a11, a21, a12, a22, alpha, psi, X, theta, knot, tvec, G, ff_Z, sigma, cutoff){
  means_regime1 <- get_mmean(alpha, X, theta, knot, tvec, a11, a21, G, ff_Z, sigma, cutoff)
  means_regime2 <- get_mmean(alpha, X, theta, knot, tvec, a12, a22, G, ff_Z, sigma, cutoff)
  var_regime1 <- 
    mapply(a1 = means_regime1[,'a1'], a2 = means_regime1[,'a2'],
         tval = means_regime1[,'time'],
         FUN = function(a1, a2, tval){
           return(get_margvar(tval=tval, alpha, psi, knot, a1=a1, a2=a2, G, ff_Z, sigma, cutoff))
         })
  var_regime2 <- 
    mapply(a1 = means_regime2[,'a1'], a2 = means_regime2[,'a2'],
           tval = means_regime2[,'time'],
           FUN = function(a1, a2, tval){
             return(get_margvar(tval=tval, alpha, psi, knot, a1=a1, a2=a2, G, ff_Z, sigma, cutoff))
           })
  
    return(
      cbind(time = means_regime1[,'time'],
          a11=means_regime1[,'a1'], a21 = means_regime1[,'a2'],
          a21 = means_regime2[,'a1'], a22=means_regime2[,'a2'],
          meandiff = means_regime1[,'mu'] - means_regime2[,'mu'],
          avgvar = (1/2) * (var_regime1 + var_regime2),
          sdavg = sqrt( (1/2) * (var_regime1 + var_regime2) ),
          effsize = ( means_regime1[,'mu'] - means_regime2[,'mu'] ) / 
            sqrt( (1/2) * (var_regime1 + var_regime2) ) )
    )
}

# theta: coefficient for X (baseline covariate)
# psi: vector of 2 coefficients, for main effect of R(a1) in second stage
#    (psi1 * indic(a1=1) + psi(-1) * indic(a1=-1)) * (R(a1) - Pr(R(a1)=1))
datfunc_mm <- function(N, G, tvec, knot, sigma, X, alpha, theta, psi, cutoff,
                               ff_Z, pr_a1 = 0.5, pr_a2 = 0.5, return_po = FALSE) {
  ni <- length(tvec)
  Nni <- N * ni
  IDs <- rep(1:N, each=ni)
  if (length(X) != N){
    stop("provide N-length X vector")
  }
  if (!(knot %in% tvec)){
    stop("knot value must be an observation time")
  }
  dim_raneff <- ncol(G)
  pir1 <- c('1' = 1 - get_pinr(alpha, knot, 1,G, ff_Z, sigma,cutoff),
            '-1' = 1 - get_pinr(alpha, knot, -1, G, ff_Z, sigma,cutoff))
  Xrep <- X[IDs]
  tlong <- rep(tvec, N)
  regimemat <- as.matrix(expand.grid(id=IDs, a1=c(1,-1), a2=c(1,-1)))
  Zlong <- model.matrix(ff_Z, data.frame(cbind(time=tlong,regimemat,Y=-99)))
  f <- function(){
    #####
    ## Generate the potential outcomes
    #####
    true_resid <- rnorm(n=Nni, mean=0, sd=sigma)
    gamma_n <- mvrnorm(n=N, mu=rep(0, dim_raneff), Sigma=G)
    gamma_rep <- gamma_n[regimemat[,'id'],,drop=FALSE]
    ranef_realized <- 
      vector('numeric', length=nrow(regimemat))
    for (j in 1:ncol(Zlong)){
      ranef_realized <- ranef_realized + Zlong[,j] * gamma_rep[,j]
    }
    
    dat <- 
      cbind(
        time=tlong,
        true_resid = true_resid, 
        regimemat,
        X=Xrep,
        Ya1a2 = as.numeric(rep(NA, Nni)),
        Ra1 = as.numeric(rep(NA,  Nni)),
        ranef = ranef_realized)
    
    dat_stage1 <- dat[dat[,'time'] <= knot,]
    
    Y_stage1 <- alpha[1] + alpha[2] * dat_stage1[,'time'] + 
      alpha[3] * dat_stage1[,'a1'] * dat_stage1[,'time'] + 
      dat_stage1[,'ranef'] + 
      theta * dat_stage1[,'X'] + 
      dat_stage1[,'true_resid'] 
    
    dat[dat[,'time'] <= knot, 'Ya1a2'] <- Y_stage1
    
    Ra1 <- 1 * ( dat[dat[,'time'] == knot, 'Ya1a2'] - theta * dat[dat[,'time']==knot,'X'] > cutoff )
    dat[,'Ra1'] <- rep(Ra1, each=ni) 
    
    dat_stage2 <- dat[dat[,'time'] > knot,]
    
    Y_stage2 <- 
      alpha[1] + alpha[2] * knot + alpha[3] * knot * dat_stage2[,'a1'] +
      (dat_stage2[,'time'] - knot) * alpha[4] + 
      (dat_stage2[,'time'] - knot) * alpha[5] * dat_stage2[,'a1'] +
      (dat_stage2[,'time'] - knot) * alpha[6] * dat_stage2[,'a2'] * (1 - dat_stage2[,'Ra1']) + 
      (dat_stage2[,'time'] - knot) * alpha[7] * dat_stage2[,'a1'] * dat_stage2[,'a2'] * (1 - dat_stage2[,'Ra1']) +
      (dat_stage2[,'time'] - knot) * (psi['1'] * 1 * (dat_stage2[,'a1']==1) + 
                                        psi['-1'] * 1 * (dat_stage2[,'a1']==-1)) * (dat_stage2[,'Ra1'] - 1 * (dat_stage2[,'a1']==1) * pir1['1'] - 1 * (dat_stage2[,'a1']==-1) * pir1['-1']) +  
      theta * dat_stage2[,'X'] + 
      dat_stage2[,'ranef'] + 
      dat_stage2[,'true_resid']
    
    dat[dat[,'time']>knot,'Ya1a2'] <- Y_stage2
    ##################
    dat_po <- dat
    
    ######################
    ## Generate the observed data using randomization and consistency
    A1 <- 2 * rbinom(n=N, size=1, prob = pr_a1) - 1
    dat <- cbind(dat,A1=A1[dat[,'id']])
    dat_obs <- dat[dat[,'a1']==dat[,'A1'],]
    dat_obs <- dat_obs[order(dat_obs[,'id'],dat_obs[,'time']),]
    
    Ra1_uniqdat <- dat_obs[dat_obs[,'time']==0 & dat_obs[,'a2']==1,]
    
    A2 <- ifelse(Ra1_uniqdat[,'Ra1']==1, NA,
                 2 * rbinom(n=1, size=1, prob=pr_a2) -1 )
    dat_obs <- cbind(dat_obs, A2 = A2[dat_obs[,'id']],
                     A2filter = ifelse(is.na(A2[dat_obs[,'id']]), 1,
                                       A2[dat_obs[,'id']]))
    
    # A2filter: 1 for Ra1=1, A2 otherwise
    #       since Ya1a2 is equal for each value of a2 when Ra1=1,
    #       we can arbitrarily pick the vector of Ya1a2 for a2=1 as the observed Y
    dat_obs <- dat_obs[dat_obs[,'A2filter']==dat_obs[,'a2'],]
    #######################
    
    ret <- 
      cbind(
        id = dat_obs[,'id'],
        time= dat_obs[,'time'],
        X = dat_obs[,'X'],
        Y = dat_obs[,'Ya1a2'],
        R = dat_obs[,'Ra1'],
        A1 = dat_obs[,'A1'],
        A2 = dat_obs[,'A2']
            )
    if (return_po) {
      return(list(obs=ret,po=dat_po))
    } else {
      return(ret)
    }
  }
  return(f)
}

get_aug_weight <- function(dobs, pr_a1= 0.5, pr_a2 = 0.5){
  d0 <- dobs[dobs[,'R']==0,]
  d1 <- dobs[dobs[,'R']==1,]
  N0 <- length(unique(d0[,'id']))
  N1 <- length(unique(d1[,'id']))
  M <- N0 + N1 + N1
  d0 <- cbind(d0, W=(1 / (pr_a1 * pr_a2)))
  d1 <- cbind(d1, W=(1 / pr_a1))
  ni <- length(unique(d0[,'time']))
  d0 <- cbind(d0, waves = rep(1:ni, N0))
  d1a <- d1
  d1b <- d1
  d1a[,'A2'] <- 1
  d1b[,'A2'] <- -1
  d1b <- cbind(d1b, waves=rep((ni+1):(ni*2), N1))
  d1a <- cbind(d1a, waves=rep(1:ni, N1))
  daw <- rbind(d0, d1a, d1b)
  daw <- daw[order(daw[,'id'],daw[,'waves']),]
  daw <- cbind(daw, idrep = rep(1:M, each=ni))
  return(daw)
}

get_2aug <- function(dobs, pr_a1=0.5, pr_a2=0.5){
  d0 <- dobs[dobs[,'R']==0,]
  d1 <- dobs[dobs[,'R']==1,]
  N0 <- length(unique(d0[,'id']))
  N1 <- length(unique(d1[,'id']))
  M <- N0 + N1 + N1
  d0 <- cbind(d0, W=4)
  d1 <- cbind(d1, W=2)
  ni <- length(unique(d0[,'time']))
  d0 <- cbind(d0, waves = rep(1:ni, N0))
  d1a <- d1
  d1b <- d1
  d1a[,'A2'] <- 1
  d1b[,'A2'] <- -1
  d1b <- cbind(d1b, waves=rep((ni+1):(ni*2), N1))
  d1a <- cbind(d1a, waves=rep(1:ni, N1))
  daw <- rbind(d0, d1a, d1b)
  daw <- daw[order(daw[,'id'],daw[,'waves']),]
  daw <- cbind(daw, idrep = rep(1:M, each=ni))
  
  daw2 <- daw[rep(1:nrow(daw), times= daw[,'W']),]
  daw2 <- cbind(daw2, 
                incr = do.call('c',sapply(daw[,'W'],function(ww) return(seq(0, (ww-1) * (1/ww), 1/ww)))))
  daw2 <- cbind(daw2, id2rep = daw2[,'idrep'] + daw2[,'incr'])
  daw2 <- daw2[order(daw2[,'id2rep'], daw2[,'time']),]
  return(
    cbind(
      id = daw2[,'id'],
      id2rep = daw2[,'id2rep'],
      time = daw2[,'time'],
      X = daw2[,'X'],
      A1 = daw2[,'A1'],
      A2 = daw2[,'A2'],
      R = daw2[,'R'],
      Y = daw2[,'Y']
    )
  )
}
