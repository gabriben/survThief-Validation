/*  Variable naming:
    obs       = observed
    cen       = (right) censored
    N         = number of samples
    M         = number of covariates
    bg        = established risk (or protective) factors
*/

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> M;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M] Xobs;
  matrix[Ncen, M] Xcen;
  int gObs[Nobs];
  int gCen[Ncen];
  /* int Xobs[Nobs, 1]; */
  /* int Xcen[Ncen, 1]; */
  /* int Xobs[Nobs]; */
  /* int Xcen[Ncen]; */
  /* int Obs[Nobs]; // make sure {1,2} */
  /* int Cen[Ncen]; */
}

parameters {
  vector[M] beta;
  real mu;
  /* real alpha[M+1]; */
  vector[2] alpha;
}

model {
  beta ~ normal(0, 5);
  alpha ~ normal(0, 5);
  mu ~ normal(0, 5);

  
  /* yobs ~ weibull(alpha[Obs], exp(-(mu + Xobs * beta)/alpha[Obs])); */
  /* target += weibull_lccdf(ycen | alpha[Cen], */
  /*                         exp(-(mu + Xcen * beta)/alpha[Cen])); */

  /* for (i in 1:Nobs) { */
  /*   yobs[i] ~ weibull(alpha[Xobs[i] + 1], */
  /*                     exp(-(mu + Xobs[i] * beta)/ */
  /*                         alpha[Xobs[i] + 1])); */
  /* } */
  /* for (i in 1:Ncen) { */
  /*   target += weibull_lccdf(ycen[i] | alpha[Xcen[i] + 1], */
  /*                           exp(-(mu + Xcen[i] * beta)/ */
  /*                               alpha[Xcen[i] + 1])); */
  /* } */
  /* int i = 1; */
  
  /* for (i in 1:Nobs) { */
  /*   yobs[i] ~ weibull(alpha[Xobs[i, 1] + 1], */
  /*                     exp(-(mu + Xobs[i, 1] * beta)/ */
  /*                         alpha[Xobs[i, 1] + 1])); */
  /* } */
  /* for (i in 1:Ncen) { */
  /*   target += weibull_lccdf(ycen[i] | alpha[Xcen[i, 1] + 1], */
  /*                           exp(-(mu + Xcen[i, 1] * beta)/ */
  /*                               alpha[Xcen[i, 1] + 1])); */
  /* } */
  
  for (i in 1:Nobs) {
    yobs[i] ~ weibull(alpha[gObs[i]],
                      exp(-(mu + Xobs[i] * beta)/
                          alpha[gObs[i]]));
  }
  for (i in 1:Ncen) {
    target += weibull_lccdf(ycen[i] | alpha[gCen[i]],
                            exp(-(mu + Xcen[i, 1] * beta)/
                                alpha[gCen[i]]));
  }
}

generated quantities {
  real yhat_uncens[Nobs + Ncen];
  real log_lik[Nobs + Ncen];
  real lp[Nobs + Ncen];
  
  for (i in 1:Nobs) {
    lp[i] = mu + Xobs[i,] * beta;
    yhat_uncens[i] = weibull_rng(alpha[gObs[i]],
                                 exp(-(mu + Xobs[i,]
                                       * beta)/alpha[gObs[i]]));
    log_lik[i] = weibull_lpdf(yobs[i] | alpha[gObs[i]],
                              exp(-(mu + Xobs[i,]
                                    * beta)/alpha[gObs[i]]));
  }

  for (i in 1:Ncen) {
    lp[Nobs + i] = mu + Xcen[i,] * beta;
    yhat_uncens[Nobs + i] = weibull_rng(alpha[gCen[i]],
                                        exp(-(mu + Xcen[i,]
                                              * beta)/alpha[gCen[i]]));
    log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha[gCen[i]],
                                      exp(-(mu + Xcen[i, 1]
                                            * beta)/alpha[gCen[i]]));
  }
}


/* generated quantities { */
/*   real yhat_uncens[Nobs + Ncen]; */
/*   real log_lik[Nobs + Ncen]; */
/*   real lp[Nobs + Ncen]; */
  
/*   for (i in 1:Nobs) { */
/*     lp[i] = mu + Xobs[i,] * beta; */
/*     yhat_uncens[i] = weibull_rng(alpha[Xobs[i, 1] + 1], */
/*                                  exp(-(mu + Xobs[i, 1] */
/*                                        * beta)/alpha[Xobs[i, 1] + 1])); */
/*     log_lik[i] = weibull_lpdf(yobs[i] | alpha[Xobs[i, 1] + 1], */
/*                               exp(-(mu + Xobs[i,] */
/*                                     * beta)/alpha[Xobs[i, 1] + 1])); */
/*   } */

/*   for (i in 1:Ncen) { */
/*     lp[Nobs + i] = mu + Xcen[i,] * beta; */
/*     yhat_uncens[Nobs + i] = weibull_rng(alpha[Xcen[i, 1] + 1], */
/*                                         exp(-(mu + Xcen[i,] */
/*                                               * beta)/alpha[Xcen[i, 1] + 1])); */
/*     log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha[Xcen[i, 1] + 1], */
/*                                       exp(-(mu + Xcen[i, 1] */
/*                                             * beta)/alpha[Xcen[i, 1] + 1])); */
/*   } */
/* } */
