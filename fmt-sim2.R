library(tidyverse)
library(xtable)
source("fmt-funcs.R")

fnames <- list.files(pattern = 'sim2-.*\\.RData')
ssizes <- unique(as.numeric(str_match(fnames, '(-N)([0-9]+)')[,3]))
effsizes <- unique(str_extract(fnames, 'small|med|large'))
effsizevalues <- c('small' = 0.2, 'med' = 0.5, 'large' = 0.8)
conflvl <- 0.95
showmethods <- c('mm_slopes','mm_intercept',
                 'exchangeable',
                 'unstructured', 
                 'independence')
                 #'independence_t', 'exchangeable_t',
                 # 'independence_t_a1a2',
                 #'ar1','exchangeable_t_a1a2')
methodkey <-
  c('exchangeable' = 'GEE Exchangeable',
    'unstructured_a1a2' = 'GEE Unstructured across time and regime',
    'unstructured' = 'GEE Unstructured',
    'ar1' = 'GEE AR',
    'mm_slopes' = 'LMM slopes and intercepts',
    'mm_intercept' = 'LMM intercepts only',
    'trueV' = 'Plug in true marginal V',
    'independence' = 'GEE Independence',
    'exchangeable_t' = 'GEE Exchangeable (variance changes over time)')

contrast_res <- NULL
vres <- NULL
for (ncur in ssizes){
  for (effsize in effsizes){
    load(list.files(pattern = paste("sim2-", effsize, 'effect-N', ncur, '-nsim', sep='')))
    contrast_res <- rbind(contrast_res, cbind('N'= ncur, '$d$' = as.numeric(effsizevalues[effsize]), 
                                          get_contrast_results(res, simparm, effsize, conflvl, showmethods)))
    vres <- rbind(vres, c('N'=ncur, '$d$' = as.numeric(effsizevalues[effsize]),
                          get_vvnorm(res, simparm, effsize, showmethods)))
  }
}


cdat <- as_tibble(contrast_res) %>% rename(d=`$d$`) %>%
  mutate(method = rownames(contrast_res))
vdat <- as_tibble(vres) %>% rename(d = `$d$`) %>%
  gather(-N, -d, key='method', value='vvnorm')

plotdat <- cdat %>% 
  left_join(vdat, by=c('N','d','method'))

plotdat$method_short <- with(plotdat, ifelse(method=='mm_slopes', 'LMM Slopes', 
                                             ifelse(method=='mm_intercept', 'LMM Intercepts Only', 'GEE')))

filter(plotdat, method %in% c('mm_slopes','mm_intercept','unstructured','independence','exchangeable')) %>%
  arrange(N, mc_rmse) %>%
  mutate(method = methodkey[method])%>%
  dplyr::select(N, Method = method,
         #`$d$` = d,
         # `True value`,
         Bias, `Monte Carlo SD` = `mc_SD`,
         `SE estimate`, `CI Coverage` = CI, `RMSE` = mc_rmse,
         `$\\|V_{\\text{true}} - \\expected{\\hat{V}}\\| / \\|V_{\\text{true}}\\|$` = vvnorm) -> tt

tt <- 
  tt %>% group_by(N)%>%
  mutate(`RMSE Inflation` = RMSE / min(RMSE)) %>% ungroup %>%
  dplyr::select(N, Method, Bias, `Monte Carlo SD`, `SE estimate`, `CI Coverage`, `RMSE Inflation`, `$\\|V_{\\text{true}} - \\expected{\\hat{V}}\\| / \\|V_{\\text{true}}\\|$`)

texFile <- file('sim2-table.tex', open = 'wt')
writeLines("% misspecification, compare to GEE", texFile)
writeLines("% medium effect size only", texFile)
invisible(capture.output(x <- print(xtable(tt, digits = c(0,0,0,3,3,3,3,3,3)),
      include.rownames = FALSE, include.colnames = TRUE, floating = FALSE,
      sanitize.colnames.function = identity, 
      math.style.negative=TRUE,
      sanitize.text.function = identity,
      hline.after=c(0,nrow(tt)),
      booktabs = TRUE)))
writeLines(x, texFile)
close(texFile)
