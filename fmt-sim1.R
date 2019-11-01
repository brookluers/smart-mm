library(tidyverse)
library(purrr)
library(xtable)
source("fmt-funcs.R")
fnames <- list.files(pattern = 'sim1-.*\\.RData')
ssizes <- unique(as.numeric(str_match(fnames, '(-N)([0-9]+)')[,3]))
effsizes <- unique(str_extract(fnames, 'small|med|large'))
effsizevalues <- c('small' = 0.2, 'med' = 0.5, 'large' = 0.8)
conflvl <- 0.95
contrast_res <- NULL
vres <- NULL
showmethods <- c('mm_slopes','mm_intercept')
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
nsimcheck <- NULL
for (ncur in ssizes){
  for (effsize in effsizes){
    load(list.files(pattern = paste("sim1-", effsize, 'effect-N', ncur, '-nsim', sep='')))
    nsimcheck <- c(nsimcheck, simparm$nsim)
    contrast_res <- rbind(contrast_res, cbind('N'= ncur, '$d$' = as.numeric(effsizevalues[effsize]), 
                                              get_contrast_results(res, simparm, effsize, conflvl, showmethods)))
    
  }
}
if (length(unique(nsimcheck)) > 1 ){
  cat("loaded .RData files with different numbers of simulations:")
  cat(unique(nsimcheck))
}
nsim <- unique(nsimcheck)


cdat <- as_tibble(contrast_res) %>% rename(d=`$d$`) %>%
  mutate(method = rownames(contrast_res))

cdat %>% filter(d %in% c(0.2,0.8)) %>%
  arrange(desc(method), d, N, mc_rmse) %>%
  mutate(method = methodkey[method])%>%
  dplyr::select(Method = method, `$d$` = d,
         `True value`,
         N,
         Bias, `Monte Carlo SD` = `mc_SD`,
         `SE estimate`, `CI Coverage` = CI, `RMSE` = mc_rmse) -> tt

texFile <- file('sim1-table.tex', open = 'wt')
writeLines("% correct model, LMM slopes and intercepts only", texFile)
invisible(capture.output(x <- print(xtable(tt, digits = c(0,0,1,3,0,3,3,3,3,3)),
                                    include.rownames = FALSE, include.colnames = TRUE, floating = FALSE,
                                    sanitize.colnames.function = identity, 
                                    math.style.negative=TRUE,
                                    sanitize.text.function = identity,
                                    hline.after=c(0,nrow(tt)),
                                    booktabs = TRUE)))

writeLines(x, texFile)
writeLines(paste("\\caption{Estimation of end-of-study contrast with slopes-and-intercepts mixed model. Values computed with ",
                 format(nsim, big.mark=',')," simulation replicates.}", sep=''), texFile)
close(texFile)
