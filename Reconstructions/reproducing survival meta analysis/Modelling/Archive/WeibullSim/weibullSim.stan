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
  real mu;
  real alpha;
}

model {
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
  mu ~ normal(0, 10);
  yobs ~ weibull(alpha, exp(-(mu + Xobs * beta)/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-(mu + Xcen * beta)/alpha));
}

generated quantities {
  real yhat_uncens[Nobs + Ncen];
  real log_lik[Nobs + Ncen];
  real lp[Nobs + Ncen];
  
  for (i in 1:Nobs) {
    lp[i] = mu + Xobs[i,] * beta;
    yhat_uncens[i] = weibull_rng(alpha, exp(-(mu + Xobs[i,]
                                              * beta)/alpha));
    log_lik[i] = weibull_lpdf(yobs[i] | alpha, exp(-(mu + Xobs[i,]
                                                     * beta)/alpha));
  }

  for (i in 1:Ncen) {
    lp[Nobs + i] = mu + Xcen[i,] * beta;
    yhat_uncens[Nobs + i] = weibull_rng(alpha, exp(-(mu + Xcen[i,]
                                                     * beta)/alpha));
    log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha, exp(-(mu + Xcen[i,]
                                                             * beta)/alpha));
  }
}
