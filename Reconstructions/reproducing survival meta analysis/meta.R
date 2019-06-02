setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Reconstructed")

studies <- list.dirs(getwd(), recursive = F, full.names = F)

## each dataset is saved in an object "img" in .RData format.
## loading RData does not allow for reassigning objects to other names, so img is kept all along.
load(paste0(studies[1], "/", studies[1], ".RData"))
d <- cbind(img$survData, study = studies[1])
for(i in studies[-1]){
  load(paste0(i, "/", i, ".RData"))
  d <- rbind(d, cbind(img$survData, study = i))
}
d <- d[order(d$time),]

setwd("../")
## save(d, file = "meta.RData")
load("meta.RData")

## KM

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=d)

library(survminer)
## pdf("results/detectCensoring V2.pdf",width=14,height=10)
ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
           risk.table.col = "strata", RUEpval = TRUE, 
           ggtheme = theme_bw(), risk.table.y.text=F,
           break.time.by=12)


#### Weibull Bayesian
## https://rstudio-pubs-static.s3.amazonaws.com/435225_07b4ab5afa824342a4680c9fb2de6098.html

## clang -> on terminal : xcode-select --install
## https://mac-how-to.gadgethacks.com/how-to/install-command-line-developer-tools-without-xcode-0168115/

library(rstan)


## simulate #################

simulWeib <- function(N, lambda, rho, beta, rateC)
{
  ## covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))

  ## Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(x * beta)))^(1 / rho)

  ## censoring times
  C <- rexp(n=N, rate=rateC)

  ## follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  ## data set

  data.frame(id=1:N,
             time=time,
             status=status,
             x=x)
}

d <- simulWeib(N=100, lambda=0.01, rho=1, beta=-0.6, rateC=0.01)

library(simsurv)

## shape: gamma <- alpha

x <- sample(x=c(0, 1), size=100, replace=TRUE, prob=c(0.5, 0.5))


covs <- data.frame(id = 1:100, trt = stats::rbinom(100, 1L, 0.5))
s1 <- simsurv(lambdas = 1, gammas = 2, betas = c(trt = -0.5),
              x = covs, maxt = 5)
d <- simsurv(lambdas = 1, gammas = 2, x = as.data.frame(sample(0:1, 100, replace = T)))



stanData <-
  list(
    ## Number of event individuals
    Nobs = sum(s1$status == 1),
    ## Number of censored individuals
    Ncen = sum(s1$status == 0),
    ## Number of covariates
    M_bg = 1,
    ## Times for event individuals
    yobs = s1$time[s1$status == 1],
    ## Times for censored individuals
    ycen = s1$time[s1$status == 0],
    ## Covariates for event individuals as a matrix
    Xobs_bg = matrix(as.numeric(s1$x == 1)[s1$status == 1]),
    ## Covariates for censored individuals as a matrix
    Xcen_bg = matrix(as.numeric(s1$x == 1)[s1$status == 0])
  )



#################

stanData <-
  list(
    ## Number of event individuals
    Nobs = sum(d$status == 1),
    ## Number of censored individuals
    Ncen = sum(d$status == 0),
    ## Number of covariates
    M = 1,
    ## Times for event individuals
    yobs = d$time[d$status == 1],
    ## Times for censored individuals
    ycen = d$time[d$status == 0],
    ## Covariates for event individuals as a matrix
    Xobs = matrix(as.numeric(d$stratum == "open")[d$status == 1]),
    ## Covariates for censored individuals as a matrix
    Xcen = matrix(as.numeric(d$stratum == "open")[d$status == 0])
  )
stanData


############


nchains <- 2
niter <- 1000
nwarmup <- 500

setwd("../")
stanWeibFit <- stan(file = "weibull.stan", data = stanData,
                    chains = nchains, iter = niter, warmup = nwarmup,
                    control=list(adapt_delta=0.99, max_treedepth = 15))

## Similarly, we can tell Stan to take smaller steps around the posterior distribution, which (in some but not all cases) can help avoid these divergences.
stanWeibFit <- stan(file = "weibull.stan", data = stanData, control = list(adapt_delta = 0.999))

save(stanWeibFit, file = "weibFit.RData")
load("weibFit.RData")

library(shinystan)
launch_shinystan(as.shinystan(stanWeibFit))

traceplot(stanWeibFit, par = c("alpha","mu","beta_bg"))

library(bayesplot)
mcmc_acf(as.matrix(stanWeibFit), pars = c("alpha","mu","beta_bg[1]"))
mcmc_areas(as.matrix(stanWeibFit), pars = c("alpha","mu","beta_bg[1]"), prob = 0.95)

library(tidybayes)
stanWeibFitDraws <- tidybayes::tidy_draws(stanWeibFit)


## divergence diagnostics
p <- c('tau_s', 'tau', 'mu', 'alpha', 'beta')
draws <- as.array(stanWeibFit, pars = p)
mcmc_parcoord(draws, transform = function(x) {(x - mean(x)) / sd(x)}, np = nuts_params(stanWeibFit))
pairs(stanWeibFit, pars = p)
print(broom::tidy(stanWeibFit, head=10))


## stanWeibFitDraws <- stanWeibFitDraws[sample(nrow(stanWeibFitDraws), 100),]

## Constructor for treatment-specific survival function
construct_survival_function <- function(alpha, mu, beta, x) {
  function(t) {
    sigma_i <- exp(-1 * (mu + beta * x) / alpha)
    exp(- (t / sigma_i)^alpha)
  }
}

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

times <- seq(from = 0, to = 160, by = 0.1)
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

## Average on survival scale
stanWeibFitSurvMean <-
  stanWeibFitSurv %>%
  group_by(treatment, t) %>%
  summarize(survival_mean = mean(survival),
            survival_95upper = quantile(survival, probs = 0.975),
            survival_95lower = quantile(survival, probs = 0.025))

save(stanWeibFitSurvMean, file = "weibFit.rds")
load("weibFit.rds")

## pdf("weibFit.pdf",width=14,height=10)
library(tidyverse)
library(ggplot2)
bayesPlot <- ggplot(data = stanWeibFitSurvMean,
                    mapping = aes(x = t, group = treatment, color = fct_rev(treatment))) +
  geom_line(aes(y = survival_mean)) +
  geom_line(aes(y = survival_95upper), linetype = "dotted") +
  geom_line(aes(y = survival_95lower), linetype = "dotted") +
  ## facet_grid(. ~ treatment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())

## add other plot

colnames(d)[1] <- "time"

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=d)

library(survminer)
## pdf("results/detectCensoring V2.pdf",width=14,height=10)
KMPlot <- ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
                     risk.table.col = "strata", RUEpval = TRUE, 
                     ggtheme = theme_bw(), risk.table.y.text=F,
                     break.time.by=12)

KMPlot$plot + geom_line(data = stanWeibFitSurvMean,
                        mapping = aes(x = t,y = survival_mean, group = treatment, color = treatment)) + 
  geom_line(data = stanWeibFitSurvMean,
            mapping = aes(x = t,y = survival_95upper, group = treatment, color = treatment), linetype = "dotted") + 
  geom_line(data = stanWeibFitSurvMean,
            mapping = aes(x = t,y = survival_95lower, group = treatment, color = treatment), linetype = "dotted")



## plot all realisations (not recommended)
ggplot(data = stanWeibFitSurv,
       mapping = aes(x = t, y = survival, color = fct_rev(treatment), group = interaction(.chain,.draw,treatment))) +
  geom_line(size = 0.1, alpha = 0.02) +
  geom_line(data = stanWeibFitSurvMean,
            mapping = aes(y = survival_mean, group = treatment)) +
  geom_line(data = stanWeibFitSurvMean,
            mapping = aes(y = survival_95upper, group = treatment),
            linetype = "dotted") +
  geom_line(data = stanWeibFitSurvMean,
            mapping = aes(y = survival_95lower, group = treatment),
            linetype = "dotted") +
  ## facet_grid(. ~ treatment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())
                                        # dev.off()


## DPWeibull


library(DPWeibull)

y <- as.matrix(d[1:2])
weibModel <- dpweib(y ~ 1)

library(fastDummies)

x <- as.matrix(cbind(1*(d$stratum == "open"), dummy_cols(d$study, remove_first_dummy = T)[-1]))

                                        # data <- data.frame(y = y, stratum = 1*(d$stratum == "open"), )

weibModelCovs <- dpweib(y ~ x)

## save(weibModelCovs, file = "../Reconstructions/reproducing survival meta analysis/DDP.RData")

load("DDP.RData")

str(weibModelCovs)

myWeib <- function(y, a, l) return(a * l * y^(a-1) * exp(-l * y^a))

ec <- 2
G_0 <- function(n) rnorm(n, 0, 10)
n <- 100
b <- rbeta(n, 1, c)
p <- numeric(n)
p[1] <- b[1]
p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
y <- G_0(n)
theta <- sample(y, prob = p, replace = TRUE)

plot(density(theta))


## Stan DP


library(rstan)

stanData <-
  list(
    ## Number of event individuals
    Nobs = sum(d$status == 1),
    ## Number of censored individuals
    Ncen = sum(d$status == 0),
    ## Number of covariates
    M_bg = 1,
    ## Times for event individuals
    yobs = d$time[d$status == 1],
    ## Times for censored individuals
    ycen = d$time[d$status == 0],
    ## Covariates for event individuals as a matrix
    Xobs_bg = matrix(as.numeric(d$stratum == "open")[d$status == 1]),
    ## Covariates for censored individuals as a matrix
    Xcen_bg = matrix(as.numeric(d$stratum == "hybrid")[d$status == 0]),
    ## number of clusters
    C <- 10
  )
stanData

setwd("../")
stanWeibFit <- stan(file = "DDP.stan", data = stanData)

save(stanWeibFit, file = "DDP.RData")
load("DDP.RData")
