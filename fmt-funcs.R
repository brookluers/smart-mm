library(purrr)
get_single_coef_results <- function(coefname, res, simparm, effsize, conflvl, showmethods) {
  nregime <- length(simparm$regimenames)
  tmax <- max(simparm$tvec)
  
  coef_trueval <- simparm$truecoefs[[effsize]][[coefname]]
  coef_ix <- which(names(simparm$truecoefs[[effsize]]) == coefname)
  methodnames <- rownames(res[[1]]$bhat)
  nsim <- length(res)
  res_1coef <- array(0, dim = c(length(methodnames), 3, nsim))
  zmult <- qnorm((1-conflvl)/2,lower.tail=F)
  for (i in seq_along(res)){
    est <- res[[i]]$bhat[, coef_ix]
    se <- sapply(res[[i]]$vlist, function(vb) return(sqrt(diag(vb)[coef_ix])))
    ci_cover <- (coef_trueval <= est + zmult * se) & (coef_trueval >= est - zmult * se)
    res_1coef[ , , i] <- cbind(est, se, ci_cover)
  }
  mc_sd <- apply(res_1coef, c(1,2), sd, na.rm = TRUE)[,1]
  mc_mean <- apply(res_1coef, c(1,2), mean, na.rm = TRUE)
  mc_rmse <- apply(res_1coef, c(1,2), function(x) return(sqrt(mean((x - coef_trueval)^2))))[,1]
  res_table <- cbind(mc_mean, mc_sd, mc_rmse)
  colnames(res_table) <- c("Estimate", "SE estimate", "CI", "mc_SD", "mc_rmse")
  rownames(res_table) <- methodnames
  return(cbind('True value' = coef_trueval, 
               res_table[showmethods, c('Estimate', 'mc_SD', 'SE estimate', 'CI', 'mc_rmse'), drop=F]))
}

get_vvnorm <- function(res, simparm, effsize, showmethods){
  nregime <- length(simparm$regimenames)
  tmax <- max(simparm$tvec)
  methodnames <- rownames(res[[1]]$bhat)
  nsim <- length(res)
  
  avg_vtrue <- (1/nregime) * reduce(simparm$Vtruelist[[effsize]], `+`)
  res_vhat <- array(0, dim = c(dim(res[[1]]$vhatmat$`1.1`), nsim))
  for (i in seq_along(res)){
    res_vhat[,,i] <- (1/nregime) * reduce(res[[i]]$vhatmat,`+`)
  }
  vtrue_ltri <- avg_vtrue[lower.tri(avg_vtrue, diag=T)]
  vhat_mean <- apply(res_vhat, c(1,2), mean)
  ret <- apply(vhat_mean, 1, function(vhj) return(norm(cbind(vhj - vtrue_ltri), type='F')))
  names(ret) <- methodnames
  return(ret[showmethods])
}
  
get_contrast_results <- function(res, simparm, effsize, conflvl, showmethods){
  nregime <- length(simparm$regimenames)
  tmax <- max(simparm$tvec)
  Lmat1 <-
    model.matrix(
      simparm$ff_fixef,
      data.frame(
        Y = -99,
        time = tmax,
        A1 = 1,
        A2 = -1,
        knot = simparm$knot,
        X = simparm$X[1]
      )
    )
  Lmat2 <-
    model.matrix(
      simparm$ff_fixef,
      data.frame(
        Y = -99,
        time = tmax,
        A1 = -1,
        A2 = -1,
        knot = simparm$knot,
        X = simparm$X[1]
      )
    )
  L_contrast <- Lmat1 - Lmat2
  contrast_trueval <- as.numeric(L_contrast %*% simparm$truecoefs[[effsize]])
  methodnames <- rownames(res[[1]]$bhat)
  nsim <- length(res)
  res_contrast <- array(0, dim = c(length(methodnames), 4, nsim))
  zmult <- qnorm((1-conflvl)/2,lower.tail=F)
  for (i in seq_along(res)){
    est <- tcrossprod(L_contrast, res[[i]]$bhat)[1,]
    bias <- est - contrast_trueval
    se <- sapply(res[[i]]$vlist, function(vb) return(sqrt(L_contrast %*% vb %*% t(L_contrast))))
    ci_cover <- (contrast_trueval <= est + zmult * se) & (contrast_trueval >= est - zmult * se)
    res_contrast[ , , i] <- cbind(bias, est, se, ci_cover) # cbind(est, se, ci_cover)
  }
  mc_sd <- apply(res_contrast, c(1,2), sd, na.rm = TRUE)[,1]
  mc_mean <- apply(res_contrast, c(1,2), mean, na.rm = TRUE)[,-2]
  mc_rmse <- apply(res_contrast, c(1,2), function(x) return(sqrt(mean((x - contrast_trueval)^2))))[,2]
  res_table <- cbind(mc_mean, mc_sd, mc_rmse)
  colnames(res_table) <- c('Bias','SE estimate', 'CI', 'mc_SD', 'mc_rmse') #c("Estimate", "SE estimate", "CI", "mc_SD", "mc_rmse")
  rownames(res_table) <- methodnames
  return(cbind('True value' = contrast_trueval, res_table[showmethods, c('Bias', 'mc_SD', 'SE estimate', 'CI', 'mc_rmse'), drop=F]))
}

get_vvnorm <- function(res, simparm, effsize, showmethods){
  nregime <- length(simparm$regimenames)
  tmax <- max(simparm$tvec)
  methodnames <- rownames(res[[1]]$bhat)
  nsim <- length(res)
  
  avg_vtrue <- (1/nregime) * reduce(simparm$Vtruelist[[effsize]], `+`)
  res_vhat <- array(0, dim = c(dim(res[[1]]$vhatmat$`1.1`), nsim))
  for (i in seq_along(res)){
    res_vhat[,,i] <- (1/nregime) * reduce(res[[i]]$vhatmat,`+`)
  }
  vtrue_ltri <- avg_vtrue[lower.tri(avg_vtrue, diag=T)]
  vhat_mean <- apply(res_vhat, c(1,2), mean)
  vhat_mean_full <- 
    lapply(1:nrow(vhat_mean), function(ix) {
      vhj <- vhat_mean[ix,]
    M <- matrix(0, nrow=length(simparm$tvec), ncol=length(simparm$tvec))
    M[lower.tri(M, diag=TRUE)] <- vhj
    Mdiag <- diag(M)
    diag(M) <- 0
    M <- M + t(M)
    diag(M) <- Mdiag
    return(M)
  })
  ret <- sapply(vhat_mean_full, function(vhj) return(norm(vhj - avg_vtrue, type='F') / norm(avg_vtrue, type='F')))
  # ret <- apply(vhat_mean, 1, function(vhj) return(norm(cbind(vhj - vtrue_ltri), type='F')))
  names(ret) <- methodnames
  return(ret[showmethods])
}
