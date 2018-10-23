load('res-fmt.RData')

sumtable_byN <-
  lapply(res_byN,
       function(xN) {
         xN %>%
           group_by(method, coef, effsize) %>%
           summarise(
             mest = mean(est, na.rm = T),
             truth = first(trueval),
             sdest = sd(est, na.rm = TRUE),
             rmse = sqrt(mean((est - trueval) ^ 2, na.rm = TRUE)),
             cprob = sum(cover, na.rm = TRUE) / (n() - sum(is.nan(est))),
             vnormsq = mean(vdiffnorm ^ 2, na.rm = TRUE),
             nna = sum(is.nan(est)),
             n = n()
           ) %>% ungroup
       })

View(subset(sumtable_byN[['5000']],coef%in%c('beta2','beta5')))
