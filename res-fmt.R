library(tidyverse)
library(stringr)

simnames <- c('sim1','sim2')
effsizes <- c('small','med','large')
ssizes <- c(50,200,1000,5000) # unique(as.numeric(str_match(allfiles, '(N)([0-9]+)')[,3]))
conflvl <- 0.95
zmult <- qnorm((1-conflvl) / 2,lower.tail=F)
res_objnames <- setNames(vector('character', length(simnames)), simnames)
allfiles <- NULL


sim_effsizes <-
  setNames(lapply(simnames, function(ss)
    return(setNames(
      vector('list', length(effsizes)), effsizes
    ))),
    simnames)
sim_params <-
  setNames(lapply(simnames,
                  function(ss)
                    return(setNames(
                      lapply(effsizes,
                             function(ee)
                               return(setNames(
                                 vector('list', length(ssizes)), ssizes
                               ))), effsizes
                    ))),
           simnames)
saveobj <-
  c('simnames',
    'ssizes',
    'conflvl',
    'zmult',
    'allfiles',
    'sim_effsizes',
    'sim_params')

for (cursim in simnames) {
  res_byN <- vector('list', length=length(ssizes))
  names(res_byN) <- as.character(ssizes)
  cfiles <- list.files(pattern = paste(cursim, '-.+\\.RData', sep=''))
  allfiles <- c(allfiles, cfiles)
  res_objnames[cursim] <- paste('res_byN_', cursim, sep='')
  saveobj <- c(saveobj, res_objnames[cursim])
  for (ncur in ssizes) {
    nfnames <- cfiles[str_detect(cfiles, paste('N', ncur, '-', sep=''))]
    allres <-
      vector('list', length(nfnames)) ## a list of the results for each effect size, fixed sample size
    avg_trueVs <- vector('list', length(nfnames))
    for (i in seq_along(allres)) {
      load(nfnames[i])
      
      ni <- length(simparm$tvec)
      
      vhatnames <- paste('v',
                         paste(row(matrix('',nrow=ni,ncol=ni)),
                               col(matrix('',nrow=ni,ncol=ni)), sep=''),sep='')
      cur_effsize <- str_extract(nfnames[i], 'small|med|large')
      
      avg_trueVs[[i]] <- matrix(0, nrow=ni, ncol=ni)
      nregime <- length(simparm$Vtruelist[[cur_effsize]])
      for (cregime in names(simparm$Vtruelist[[cur_effsize]])) {
        avg_trueVs[[i]] <- avg_trueVs[[i]] + (1 / nregime) * simparm$Vtruelist[[cur_effsize]][[cregime]]
      }
      
      sim_effsizes[[cursim]][[cur_effsize]] <- get(paste(cur_effsize, '_effsizes', sep=''))
      allres[[i]] <-
        list(truecoefs = get(paste("truecoefs_", cur_effsize, sep = '')),
             res = get(paste("res_", cur_effsize, sep = '')))
      names(allres)[i] <- cur_effsize
      names(avg_trueVs)[i] <- cur_effsize
      sim_params[[cursim]][[cur_effsize]][[as.character(simparm$N)]] <- simparm
    }
    tmax <-
      max(simparm$tvec) # assuming all files with different sample sizes have the same simulation params
    Lmat1 <-
      model.matrix(
        simparm$ff_fixef,
        data.frame(
          Y = -99,
          time = tmax,
          A1 = 1,
          A2 = 1,
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
          A2 = 1,
          knot = simparm$knot,
          X = simparm$X[1]
        )
      )
    
    allres_fmt <- setNames(vector('list', length=length(effsizes)), effsizes)
    for(ceffsize in effsizes){ # for each of the effect sizes at a given sample size
      x <- allres[[ceffsize]]
      true_eos_difference <-
        Lmat1 %*% x$truecoefs - Lmat2 %*% x$truecoefs
      
      allres_fmt[[ceffsize]] <- 
        bind_rows(lapply(x$res, function(sim_i) {
          ## For each of the simulated data sets
          colnames(sim_i$bhat) <- names(x$truecoefs)
          # Gather the estimates of the regression coefficients,
          # with identifying info like N, method name and ||Vhat - V ||^2
          fmt_i <-
            data.frame(
              coef = rep(names(x$truecoefs), each = nrow(sim_i$bhat)),
              est = as.numeric(sim_i$bhat),
              method = rep(sim_i$method, length(x$truecoefs)),
              N = ncur,
              vdiffnorm = rep(sim_i$vnorms, length(x$truecoefs)),
              trueval = rep(x$truecoefs, each = length(sim_i$method))
            )
          ##
          ## Do the same for the estimate of the end-of-study contrast
          contrast_results <-
            data.frame(
              coef = 'contrast_estimate',
              vdiffnorm = sim_i$vnorms,
              est = as.numeric(Lmat1 %*% t(sim_i$bhat) - Lmat2 %*% t(sim_i$bhat)),
              N = ncur,
              method = sim_i$method,
              trueval = as.numeric(true_eos_difference)
            )
          ##
          ## Get the estimates of the marginal variance matrix
          ##
          vhat_results <-
            data.frame(
              coef = rep(vhatnames, each=nrow(sim_i$vhatmat)),
              est = as.numeric(sim_i$vhatmat),
              trueval = rep(as.numeric(avg_trueVs[[ceffsize]]), each = nrow(sim_i$vhatmat)),
              N = ncur, method = rep(sim_i$method, length(vhatnames))
            )
          
          fmt_i <- bind_rows(fmt_i, contrast_results, vhat_results)
          
          ##
          # Collate the standard errors for each coefficient and the contrast
          sedat_contrast <-
            data.frame(
              sehat = sqrt(sapply(sim_i$vlist, function(vm)
                return((Lmat1 - Lmat2) %*% vm %*% t(Lmat1 - Lmat2)
                ))),
              method = sim_i$method,
              coef = 'contrast_estimate'
            )
          sedat_beta <-
            data.frame(
              sehat = as.numeric(sapply(sim_i$vlist, function(vm)
                return(
                  setNames(sqrt(diag(vm)), names(x$truecoefs))
                ))),
              method = rep(sim_i$method, each = length(x$truecoefs)),
              coef = rep(names(x$truecoefs), length(sim_i$method))
            )
          
          sedat <- bind_rows(sedat_contrast, sedat_beta)
          
          ## Join the coefficient estimates with their standard errors
          fmt_i <-
            fmt_i %>% left_join(sedat,
                                by = c('method', 'coef')) %>%
            mutate(
              cilwr = est - zmult * sehat,
              ciupr = est + zmult * sehat,
              cover = 1 * ((cilwr <= trueval) &
                             (trueval <= ciupr))
            )
          return(fmt_i)
        }))
      
      } # end loop over effect sizes within a given sample size
    
    allres_fmt <- bind_rows(allres_fmt, .id='effsize')
    res_byN[[as.character(ncur)]] <- allres_fmt
  
    } # end loop over each sample size
  assign(res_objnames[cursim], res_byN)
}

dateFile <- file("date-formatted.tex")
writeLines(paste("Simulation results obtained on ",
      format(Sys.Date(), "%B %d, %Y"),
      " at ",
      format(Sys.time(), "%H:%M"),
      "\\\\", sep=''), dateFile)
close(dateFile)
save(list=saveobj, file = 'res-fmt.RData')
