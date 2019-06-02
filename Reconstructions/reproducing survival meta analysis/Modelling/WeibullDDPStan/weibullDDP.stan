/*  Variable naming:
    obs       = observed
    cen       = (right) censored
    N         = number of samples
    M         = number of covariates
    bg        = established risk (or protective) factors
    tau       = scale parameter
    C         = number of mixture clusters
*/

data {
  int<lower=0> C;
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> M;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M] Xobs;
  matrix[Ncen, M] Xcen;
}

parameters {
  vector[M] beta[C];
  real mu[C];
  real<lower=0> alpha[C];
  /* real mu_cl[C]; //cluster mean */
  real <lower=0,upper=1> v[C];
  /* real ps[C]; */
}

/* https://ecosang.github.io/blog/study/dirichlet-process-with-stan/ */
transformed parameters{
  simplex [C] pi;
  pi[1] = v[1];
  // stick-break process based on The BUGS book Chapter 11 (p.294)
  for(j in 2:(C-1)){
    pi[j]= v[j]*(1-v[j-1])*pi[j-1]/v[j-1]; 
  }
  pi[C]=1-sum(pi[1:(C-1)]); // to make a simplex.
}

model {
  real ps[C];
  vector[C] log_pi = log(pi);
  /* to_vector(beta) ~ normal(0,10); */
  for(c in 1:C){
    beta[c] ~ normal(0, 10);
    }
  alpha ~ normal(0, 10);
  mu ~ normal(0, 10);
  for(i in 1:Nobs){
    for(c in 1:C){
      ps[c] = log_pi[c] +
        weibull_lpdf(yobs[i] | alpha[c],
                     exp(-(mu[c] + Xobs * beta[c])/alpha[c]));
    }
    target += log_sum_exp(ps);
  }
  /* for(c in 1:C){ */
  /*   target += log(pi[c]) + weibull_lccdf(ycen[c] | alpha[c], */
  /*                                        exp(-(mu[c] + Xcen * beta[c])/alpha[c])); */
  /* } */
}

/* generated quantities { */
/*   real yhat_uncens[Nobs + Ncen]; */
/*   real log_lik[Nobs + Ncen]; */
/*   real lp[Nobs + Ncen]; */
  
/*   for (i in 1:Nobs) { */
/*     lp[i] = mu + Xobs[i,] * beta; */
/*     yhat_uncens[i] = weibull_rng(alpha, exp(-(mu + Xobs[i,] */
/*                                               * beta)/alpha)); */
/*     log_lik[i] = weibull_lpdf(yobs[i] | alpha, exp(-(mu + Xobs[i,] */
/*                                                      * beta)/alpha)); */
/*   } */

/*   for (i in 1:Ncen) { */
/*     lp[Nobs + i] = mu + Xcen[i,] * beta; */
/*     yhat_uncens[Nobs + i] = weibull_rng(alpha, exp(-(mu + Xcen[i,] */
/*                                                      * beta)/alpha)); */
/*     log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha, exp(-(mu + Xcen[i,] */
/*                                                              * beta)/alpha)); */
/*   } */
/* } */
