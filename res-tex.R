library(tidyverse)
library(stargazer)
load('res-fmt.RData')

fmtnumbers <- function(x){
  return(prettyNum(x, digits=3, big.mark=',', drop0trailing=T))
}

for (cursim in simnames){
  sumtable_byN <-
    lapply(get(paste('res_byN_', cursim, sep='')),
           function(xN) {
             xN %>%
               group_by(method, coef, effsize) %>%
               summarise(
                 mest = mean(est, na.rm = T),
                 truth = first(trueval),
                 bias = mean(est - first(trueval), na.rm=T),
                 sdest = sd(est, na.rm = TRUE),
                 rmse = sqrt(mean((est - trueval) ^ 2, na.rm = TRUE)),
                 cprob = sum(cover, na.rm = TRUE) / (n() - sum(is.nan(est))),
                 vnormsq = mean(vdiffnorm ^ 2, na.rm = TRUE),
                 nna = sum(is.nan(est)),
                 nsim = n()
               ) %>% ungroup
           })
  
  assign(paste(cursim, '_summary', sep=''),
         bind_rows(sumtable_byN, .id='n') %>% mutate(n=as.numeric(n)))
}

dispcoefs <- c('contrast_estimate', 'beta2', 'beta5','beta6')
methodkey <-
  c('exch_lucy' = 'geeglm exchangeable',
    'exch_plugin' = 'Plug-in exchangeable',
    'indep_lucy' = 'geeglm independent',
    'mm' = 'LMM Slopes and intercepts',
    'unstr_plugin' = 'Plug-in unstructured')
rawFile <- file('rawtables.txt', open='wt')  
writeLines(format(Sys.time()), con=rawFile)
for (cursim in simnames){
  writeLines(paste("\n\n---    Results for ", cursim, "    ---", sep=''), con=rawFile)
  summary_bymethod <- setNames(vector('list', length=length(methodkey)),
                               names(methodkey))
  simsummary <- get(paste(cursim, '_summary', sep=''))
  nsim <- unique(simsummary$nsim)
  writeLines(paste('\n    nsim = ', nsim, sep=''), con=rawFile)
  simsummary %>%
    mutate(n = as.numeric(n),
           effsize = factor(effsize, levels=c('small','med','large'),
                            labels=c('Effect size: small','medium','large'))) %>%
    dplyr::select(-nna, -nsim, -truth) %>%
    filter(coef %in% dispcoefs) %>%
    gather(mest, bias, sdest, rmse, cprob, vnormsq, key='measure', value='val') %>%
    mutate(val = fmtnumbers(val), method_eff=interaction(effsize, method)) %>% dplyr::select(-method,-effsize)%>%
    spread(method_eff, val) %>%
    arrange(n,coef,measure) -> simtable
    invisible(capture.output(x <- stargazer(simtable, 
                                            title=cursim,
                                            summary=F,
                                            type='text')))
    writeLines(x, con=rawFile)
}
close(rawFile)


filter(sim1_summary, coef =='beta2',method!='unstr_plugin') %>%
  mutate(n = as.numeric(n))%>%
  #filter(n>1000)%>%
  ggplot(aes(x=vnormsq, y=rmse, color=as.factor(n))) + 
#  geom_text(aes(label=method)) +
  geom_point(aes(shape=method))+
  facet_wrap(~effsize+n,scales='free')

filter(sim2_summary, coef =='beta2',method!='unstr_plugin') %>%
  mutate(n = as.numeric(n))%>%
  #filter(n>1000)%>%
  ggplot(aes(x=vnormsq, y=rmse, color=as.factor(n))) + 
  #  geom_text(aes(label=method)) +
  geom_point(aes(shape=method))+
  facet_wrap(~effsize+n,scales='free')
