/*  Variable naming:
    obs       = observed
    cen       = (right) censored
    N         = number of samples
    M         = number of covariates
    bg        = established risk (or protective) factors
    tau       = scale parameter
*/

functions{
  // https://link.springer.com/article/10.3758/s13428-016-0746-9
  vector Ga_log(vector x, int a, int b){
    vector[num_elements(x)] prob;
    real lprob;
    prob = x^(a-1) * exp(-b*x) * b^a / tgamma(a);
    lprob = sum(log(prob));
    return lprob;
  }

  int f(int λ){
    return max(0, log(log(20) /  λ)/log(25));
      
  }
  
}

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> M_bg;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  matrix[Nobs, M_bg] Xobs_bg;
  matrix[Ncen, M_bg] Xcen_bg;
}

/* transformed data { */
/*   real<lower=0> tau_mu; */
/*   real<lower=0> tau_al; */
/*   tau_mu = 10.0; */
/*   tau_al = 10.0; */
/* } */

  parameters {
    real<lower=0> tau_s_bg_raw;
    vector<lower=0>[M_bg] tau_bg_raw;
    real alpha_raw;
    vector[M_bg] beta_bg_raw;
    real mu;

    int<lower=0> C; //num of clusters

  }

  transformed parameters {
    vector[M_bg] beta_bg;
    real alpha;
    beta_bg = bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
    alpha = exp(10 * alpha_raw);

    real alpha00; //Parameter for the base distribution of λ
    alpha00 = 1.354028;
    real alpha0; //Parameter for the base distribution of λ
    alpha0 = 0.03501257;
    real lambda00; //Parameter for the base distribution of λ
    lambda00 = 7.181247;
    real alphaalpha; //Parameter for the base distribution of α
    alphaalpha = 0.2;
    real alphalambda; // Parameter for the base distribution of α
    alphalambda = 0.1;
    real a; // Parameter for the gamma prior of the concentration parameter of DP
    a = 1;
    real b; // Parameter for the gamma prior of the concentration parameter of DP
    b = 1;
    real betasl; // Parameter for the base distribution of the regression coefficients β 
    betasl = 2.5;
  
    // https://ecosang.github.io/blog/study/dirichlet-process-with-stan/
    // see also https://stats.stackexchange.com/questions/95120/simulate-dirichlet-process-in-r
    simplex [C] pi;
    pi[1] = v[1];
    // stick-break process based on The BUGS book Chapter 11 (p.294)
    for(j in 2:(C-1)){
      pi[j ]= v[j]*(1-v[j-1])*pi[j-1]/v[j-1]; 
    }
    pi[C]=1-sum(pi[1:(C-1)]); // to make a simplex.
  
  }

  model {
    yobs ~ weibull(alpha, exp(-(mu + Xobs_bg * beta_bg)/alpha));
    target += weibull_lccdf(ycen | alpha, exp(-(mu + Xcen_bg * beta_bg)/alpha));

    beta_bg_raw ~ normal(0, 1);
    alpha_raw ~ normal(0, 1);

    /* mu ~ normal(0.0, tau_mu); */
    mu ~ normal(0, 10);

    nu ~ Ga(a,b);
    lambda0 ~ Ga(alpha00, lambda00);
    // bivariate prior for (α,λ) employing a product of two gammas 
    G0 = Ga( lambda ! ,lambda0) * Identity * alpha ! *
      Ga(alphalpha, alphalambda) * q(beta_bg_raw);
  
  }

  generated quantities {
    real yhat_uncens[Nobs + Ncen];
    real log_lik[Nobs + Ncen];
    real lp[Nobs + Ncen];

    for (i in 1:Nobs) {
      lp[i] = mu + Xobs_bg[i,] * beta_bg;
      yhat_uncens[i] = weibull_rng(alpha, exp(-(mu + Xobs_bg[i,]
                                                * beta_bg)/alpha));
      log_lik[i] = weibull_lpdf(yobs[i] | alpha, exp(-(mu + Xobs_bg[i,]
                                                       * beta_bg)/alpha));
    }
    for (i in 1:Ncen) {
      lp[Nobs + i] = mu + Xcen_bg[i,] * beta_bg;
      yhat_uncens[Nobs + i] = weibull_rng(alpha, exp(-(mu + Xcen_bg[i,]
                                                       * beta_bg)/alpha));
      log_lik[Nobs + i] = weibull_lccdf(ycen[i] | alpha, exp(-(mu + Xcen_bg[i,]
                                                               * beta_bg)/alpha));
    }
  }
