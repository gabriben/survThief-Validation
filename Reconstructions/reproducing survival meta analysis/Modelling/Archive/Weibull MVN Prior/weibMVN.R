library("rstan")

load("../../meta.RData")

d$sim <- rnorm(nrow(d), 0, 1)

stanData <-
  list(
    ## Number of event individuals
    Nobs = sum(d$status == 1),
    ## Number of censored individuals
    Ncen = sum(d$status == 0),
    ## Number of covariates
    M = 2,
    ## Times for event individuals
    yobs = d$time[d$status == 1],
    ## Times for censored individuals
    ycen = d$time[d$status == 0],
    ## Covariates for event individuals as a matrix
    Xobs = cbind(as.numeric(d$stratum == "open")[d$status == 1], (d$sim)[d$status == 1]),
    ## Covariates for censored individuals as a matrix
    Xcen = cbind(as.numeric(d$stratum == "open")[d$status == 0], (d$sim)[d$status == 0])
  )

setwd("../Weibull MVN Prior")
stanWeibFit <- stan(file = "weibullMVN.stan", data = stanData)

save(stanWeibFit, file = "weibFitMVN.RData")
load("weibFit.RData")

library(shinystan)
launch_shinystan(as.shinystan(stanWeibFit))

traceplot(stanWeibFit, par = c("alpha","mu","beta"))

library(bayesplot)
mcmc_acf(as.matrix(stanWeibFit), pars = c("alpha","mu","beta"))
mcmc_areas(as.matrix(stanWeibFit), pars = c("alpha","mu","beta"), prob = 0.95)

list_of_draws <- extract(stanWeibFit)
print(names(list_of_draws))

draws <- as.array(stanWeibFit, pars = names(list_of_draws))
mcmc_parcoord(draws, transform = function(x) {(x - mean(x)) / sd(x)}, np = nuts_params(stanWeibFit))
pairs(stanWeibFit, pars = c('tau_s_bg_raw', 'tau_bg_raw', 'alpha_raw', ' beta_bg_raw', 'mu'))

## Constructor for treatment-specific survival function
construct_survival_function <- function(alpha, mu, beta, x) {
  function(t) {
    sigma_i <- exp(-1 * (mu + beta * x) / alpha)
    exp(- (t / sigma_i)^alpha)
  }
}

library(tidybayes)
stanWeibFitDraws <- tidybayes::tidy_draws(stanWeibFit)

library(dplyr)
library(purrr)
## Random functions
stanWeibFitFunc <-
  stanWeibFitDraws %>%
                      select(.chain, .iteration, .draw, alpha, mu, `beta_bg[1]`) %>%
                      ## Simplify name
                      rename(beta = `beta_bg[1]`) %>%
                      ## Construct realization of random functions
                      mutate(`S(t|1)` = pmap(list(alpha, mu, beta), function(a,m,b) {construct_survival_function(a,m,b,1)}),
                             `S(t|0)` = pmap(list(alpha, mu, beta), function(a,m,b) {construct_survival_function(a,m,b,0)}))

up10 <- function(x) ceiling(x/10)*10

times <- seq(from = 0, to = up10(max(d$time)), by = 0.1)
times_df <- data_frame(t = times)

## Try first realizations
stanWeibFitFunc$`S(t|1)`[[1]](times[1:10])
stanWeibFitFunc$`S(t|0)`[[1]](times[1:10])

library(tidyr)

stanWeibFitSurv <-
  stanWeibFitFunc %>%
  mutate(times_df = list(times_df)) %>%
  mutate(times_df = pmap(list(times_df, `S(t|1)`, `S(t|0)`),
                         function(df, s1, s0) {df %>% mutate(s1 = s1(t),
                                                             s0 = s0(t))})) %>%
  select(-`S(t|1)`, -`S(t|0)`) %>%
  unnest() %>%
  gather(key = treatment, value = survival, s1, s0) %>%
  mutate(treatment = factor(treatment,
                            levels = c("s1","s0"),
                            labels = c("open","hybrid")))

## Average and quantiles in survival scale

probs = seq(0.025, 0.975, 0.025)
probs <- probs[!probs == 0.5]
quantiles <- ""
for(i in probs) { 
  nam <- paste("q", i, sep = "")
  # assign(nam, quantile(stanWeibFitSurv$survival, probs = i))
  quantiles <- paste(quantiles, nam, "<- quantile(survival, probs = ", i, "),")
}
quantiles <- substr(quantiles, 1, nchar(quantiles) - 1)

textCode <- paste("stanWeibFitSurvMean <- stanWeibFitSurv %>% group_by(treatment, t) %>% summarize(survival_mean = mean(survival),", quantiles, ")")

eval(parse(text=textCode))

save(stanWeibFitSurvMean, file = "weibFitGranular.rds")

## pdf("weibFit.pdf",width=14,height=10)
library(tidyverse)
library(ggplot2)

## combine KM and Weibull

colnames(d)[1] <- "time"

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=d)

library(survminer)
## pdf("results/detectCensoring V2.pdf",width=14,height=10)
KMPlot <- ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
                     risk.table.col = "strata", RUEpval = TRUE, 
                     ggtheme = theme_bw(), risk.table.y.text=F,
                     break.time.by=12)

stanWeibFitSurvMean %>% 
        .[grepl(as.character(rev(probs)[i]), names(.))]

i <- 10

combPlot <- KMPlot$plot + 
  geom_line(data = stanWeibFitSurvMean,
                        mapping = aes(x = t,y = survival_mean, group = treatment, color = treatment)) + 
  geom_ribbon(data = stanWeibFitSurvMean, 
               mapping = aes(x=t, group = treatment,
                             ymax=stanWeibFitSurvMean %>% 
                                          .[grepl(as.character(rev(probs)[i]), names(.))] %>%
                               unlist(use.names = FALSE),
                             ymin=stanWeibFitSurvMean %>% 
                                          .[grepl(as.character(probs[i]), names(.))] %>%
                               unlist(use.names = FALSE), 
               fill=cols[i], alpha=.5), inherit.aes = FALSE)

bayesPlot <- ggplot(data = stanWeibFitSurvMean,
                    mapping = aes(x = t, group = treatment, color = fct_rev(treatment))) +
  geom_line(aes(y = survival_mean)) +
  geom_ribbon(data = stanWeibFitSurvMean, 
              mapping = aes(x=t, 
                            ymax=stanWeibFitSurvMean %>% 
                              .[grepl(as.character(rev(probs)[i]), names(.))] %>%
                              unlist(use.names = FALSE),
                            ymin=stanWeibFitSurvMean %>% 
                              .[grepl(as.character(rev(probs)[i]), names(.))] %>%
                              unlist(use.names = FALSE), 
                            fill=cols[i], alpha=.5))


colfunc <- colorRampPalette(c("red", "pink"))
cols <- colfunc(length(probs)/2)

combPlot <- KMPlot$plot + 
  geom_line(data = stanWeibFitSurvMean,
            mapping = aes(x = t,y = survival_mean, group = treatment, color = treatment)) 

for(i in 1: length(probs)/2){
  combPlot <-  combPlot + 
    geom_ribbon(data = stanWeibFitSurvMean, 
                mapping = aes(x=t, group = treatment,
                              ymax=stanWeibFitSurvMean %>% 
                                .[grepl(as.character(rev(probs)[i]), names(.))] %>%
                                unlist(use.names = FALSE),
                              ymin=stanWeibFitSurvMean %>% 
                                .[grepl(as.character(probs[i]), names(.))] %>%
                                unlist(use.names = FALSE), 
                              fill=cols[i], alpha=1), inherit.aes = FALSE)
}
