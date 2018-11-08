library(tidyverse)
library(stargazer)
library(stringr)
library(xtable)
load('res-fmt.RData')

fmtnumbers <- function(x){
  return(prettyNum(x, digits=3, big.mark=',', drop0trailing=T))
}

for (cursim in simnames){
  sumtable_byN <-
    lapply(get(paste('res_byN_', cursim, sep='')),
           function(xN) {
            
             regcoefs <- 
               xN %>%
               filter(!str_detect(coef, 'v\\d\\d')) %>% 
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
            
                vhatnorms <- 
                  filter(xN, str_detect(coef, 'v\\d\\d')) %>% 
                group_by(method, coef, effsize) %>%
                summarise(mest = mean(est, na.rm=T),
                          truth = first(trueval),
                          nna = sum(is.nan(est))) %>%
                  group_by(method, effsize) %>%
                  summarise(v_v_norm = sqrt(sum((mest - truth)^2))) %>% ungroup
                
                return(
                  regcoefs %>% left_join(vhatnorms,
                                       by=c('method','effsize'))
                )
              
           })
  
  assign(paste(cursim, '_summary', sep=''),
         bind_rows(sumtable_byN, .id='n') %>% mutate(n=as.numeric(n)))
}

dispeffsizes <- c('small','large')
dispcoefs <- c('contrast_estimate', 'beta2','beta5','beta6')
dispmethods <- c('exch_plugin', 'unstr_plugin', 'indep_lucy', 'mm_slopes','mm_intercept','trueV')
methodkey <-
  c('exch_lucy' = 'geeglm exchangeable',
    'exch_plugin' = 'Plug-in exchangeable',
    'indep_lucy' = 'geeglm independent',
    'mm_slopes' = 'LMM Slopes and intercepts',
    'mm_intercept' = 'LMM Intercept only',
    'trueV' = 'Plug in true marginal V',
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
    filter(effsize %in% dispeffsizes) %>%
    mutate(n = as.numeric(n),
           effsize = factor(effsize, levels=c('small','med','large'),
                            labels=c('Effect size: small','medium','large'))) %>%
    dplyr::select(-nna, -nsim, -truth, -vnormsq) %>%
    filter(coef %in% dispcoefs, method %in% dispmethods) %>%
    gather(mest, bias, sdest, rmse, cprob, v_v_norm, key='measure', value='val') %>%
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


dispeffsizes <- c('small','large')
dispcoefs <- c('contrast_estimate')
dispmethods <- c('unstr_plugin', 'exch_plugin', 'indep_lucy', 'mm_slopes')
methodkey <-
  c('exch_lucy' = 'geeglm exchangeable',
    'exch_plugin' = 'Plug-in exchangeable',
    'indep_lucy' = 'geeglm independent',
    'mm_slopes' = 'LMM Slopes and intercepts',
    'mm_intercept' = 'LMM Intercept only',
    'trueV' = 'Plug in true marginal V',
    'unstr_plugin' = 'Plug-in unstructured')
measurekey <- c(
  'v_v_norm' = '$\\|V^* - V\\|$',
  'rmse' = 'RMSE',
  'cprob' = 'CI Coverage',
  'bias' = 'Bias',
  'sdest' = 'SD'
)

captions <-
  c('sim1' = '\\caption{Slopes and intercepts generative model, $\\alpha_5==\\alpha_6=\\psi^{(1)}=\\psi^{(-1)}=0$. Values are computed from ',
    'sim2' = '\\caption{Slopes and intercepts generative model, $\\alpha_5$ and other generative regression parameters are nonzero. Values are computed from ')
effsizekey <- c(
  'large' = '$d\\approx 0.8$',
  'small' = '$d\\approx 0.2$',
  'med' = '$d\\approx 0.5$'
)
for (cursim in simnames){
  texFile <- file(paste(cursim, '-table.tex', sep=''), open='wt')  
  # writeLines(format(Sys.time()), con=rawFile)
  summary_bymethod <- setNames(vector('list', length=length(methodkey)),
                               names(methodkey))
  simsummary <- get(paste(cursim, '_summary', sep=''))
  nsim <- unique(simsummary$nsim)
  captions[cursim] <- paste(captions[cursim], nsim, ' simulated data sets.}', sep='', collapse='')
  simsummary %>%
    filter(effsize %in% dispeffsizes) %>%
    mutate(n = as.numeric(n)) %>%
    dplyr::select(-nna, -nsim, -truth, -mest, -vnormsq) %>%
    filter(coef %in% dispcoefs, method %in% dispmethods) %>%
    gather(bias, sdest, rmse, cprob, v_v_norm, key='measure', value='val') %>%
    mutate(val = fmtnumbers(val), method_eff=interaction(effsize, method)) %>% 
    dplyr::select(-method,-effsize, -coef) %>%
    spread(method_eff, val) %>%
    arrange(n, measure) %>% 
    mutate(measure = measurekey[measure]) -> simtable
  
  simtable <- bind_rows(
    lapply(unique(simtable$n), function(ncur) {
    stn <- filter(simtable,n==ncur)
    bind_rows(tibble(measure=paste('$N=', ncur, '$', sep='')),
              dplyr::select(stn,-n))
  }))
  
  stfinal <- simtable[,1]
  j <- 1;
  dispmethods <- unique(do.call('rbind',str_split(colnames(simtable)[-1], '\\.'))[,2])
  for (cmethod in dispmethods){
      stfinal <- cbind(stfinal, simtable[,(j+1):(j + length(dispeffsizes))],'')
      colnames(stfinal)[ncol(stfinal)] <- paste('blank',j,sep='')
      j <- j + length(dispeffsizes)
  }
  cmrstart <- which(str_detect(colnames(stfinal), 'blank')) - length(dispeffsizes)
  stfinal <- dplyr::select(stfinal, -ncol(stfinal))
  one_method_align <- paste( paste(rep('r', length(dispeffsizes)), collapse=''), 'c', sep='')
  effsizeline <- paste(   effsizekey[do.call('rbind', str_split(colnames(stfinal)[-1], '\\.'))[,1]], collapse='&',
                          sep='')
  effsizeline <- paste('&', str_replace_all(effsizeline, 'NA', ''),
                       '\\\\\\midrule', sep='', collapse='')
  topline <- paste('\\begin{tabular}{r', 
        paste( paste(rep( one_method_align, length(dispmethods) - 1),collapse='',sep=''), 
               paste(rep('r', length(dispeffsizes)), sep='', collapse=''), sep='',collapse=''),
        '}\\toprule',
        sep='',collapse='')
  
  headerline <- 
    paste(paste('\\multicolumn{', length(dispeffsizes), '}{c}',sep='', collapse=''),
        '{', methodkey[dispmethods], '}', sep='', collapse='& \\phantom{abc}& ')
  headerline <- paste('&', headerline, '\\\\', sep='', collapse='')
  cmidrules <- 
    paste('\\cmidrule{',
    paste(cmrstart, '-', cmrstart + length(dispeffsizes)-1, sep=''),
    '}', sep='', collapse='')
  writeLines(topline, texFile)
  writeLines(headerline, texFile)
  writeLines(cmidrules, texFile)
  writeLines(effsizeline, texFile)
  invisible(capture.output(x <- print(xtable(stfinal),
        only.contents=TRUE,
        hline.after=NULL,
        sanitize.text.function=identity,
        floating=F,
        include.rownames=F, include.colnames=F)))
  writeLines(x, texFile)
  writeLines('\\bottomrule\n\\end{tabular}', texFile)
  writeLines(captions[cursim], texFile)
  close(texFile)
}


filter(sim1_summary, coef =='contrast_estimate', method != 'exch_lucy') %>%
  mutate(n = as.numeric(n))%>%
  #filter(n>1000)%>%
  ggplot(aes(x=v_v_norm, y=rmse, color=as.factor(n))) + 
#  geom_text(aes(label=method)) +
  geom_point(aes(shape=method))+
  facet_wrap(~effsize+n,scales='free')

filter(sim2_summary, coef =='contrast_estimate',method!='exch_lucy') %>%
  mutate(n = as.numeric(n))%>%
  #filter(n>1000)%>%
  ggplot(aes(x=v_v_norm, y=rmse, color=as.factor(n))) + 
  #  geom_text(aes(label=method)) +
  geom_point(aes(shape=method))+
  facet_wrap(~effsize+n,scales='free')
