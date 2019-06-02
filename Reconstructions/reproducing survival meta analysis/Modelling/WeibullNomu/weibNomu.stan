/*  Variable naming:
    obs       = observed
    cen       = (right) censored
    N         = number of samples
    M         = number of covariates
    bg        = established risk (or protective) factors
    tau       = scale parameter
*/

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> M;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M] Xobs;
  matrix[Ncen, M] Xcen;
}

parameters {
  vector[M] beta;
  real<lower=0> alpha;
}

model {
  beta ~ normal(0, 5);
  alpha ~ normal(0, 5);
  yobs ~ weibull(alpha, exp(-(Xobs * beta)/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-(Xcen * beta)/alpha));
}

generated quantities {
  real yhat_uncens[Nobs + Ncen];
  real log_lik[Nobs + Ncen];
  real lp[Nobs + Ncen];
  
  for (i in 1:Nobs) {
    lp[i] = Xobs[i,] * beta;
    yhat_uncens[i] = weibull_rng(alpha, exp(-(Xobs[i,]
                                              * beta)/alpha));
    log_lik[i] = weibull_lpdf(yobs[i] | alpha, exp(-(Xobs[i,]
                                                     * beta)/alpha));
  }

  for (i in 1:Ncen) {
    lp[Nobs + i] = Xcen[i,] * beta;
    yhat_uncens[Nobs + i] = weibull_rng(alpha, exp(-(Xcen[i,]
                                                     * beta)/alpha));
    log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha, exp(-(Xcen[i,]
                                                             * beta)/alpha));
  }
}
