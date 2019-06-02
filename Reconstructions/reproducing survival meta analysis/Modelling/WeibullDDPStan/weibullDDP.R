setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Modelling/WeibullDDPStan/")


load("meta.RData")

library(rstan)

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
    Xcen = matrix(as.numeric(d$stratum == "open")[d$status == 0]),
    C = 10
    ##  number of mixture clusters
  )

stanWeibDDPFit <- stan(file = "weibullDDP.stan", data = stanData, chains=1, init_r = 0.01, refresh =1)
