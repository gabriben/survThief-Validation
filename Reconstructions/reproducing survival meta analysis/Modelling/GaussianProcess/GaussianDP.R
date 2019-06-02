rm(list = ls())
library(mixtools)
library(ggplot2)
library(tidyverse)
library(magrittr)
## Data generation code retrieved from
## http://www.jarad.me/615/2013/11/13/fitting-a-dirichlet-process-mixture

dat_generator <- function(truth) {
  set.seed(1)
  n = 500
  
  f = function(x) {
    out = numeric(length(x))
    for (i in 1:length(truth$pi)) out = out + truth$pi[i] * dnorm(x, truth$mu[i], 
                                                                  truth$sigma[i])
    out
  }
  y = rnormmix(n, truth$pi, truth$mu, truth$sigma)
  for (i in 1:length(truth$pi)) {
    assign(paste0("y", i), rnorm(n, truth$mu[i], truth$sigma[i]))
  }
  dat <- data_frame(y = y, y1 = y1, y2 = y2, y3 = y3)
}
truth = data.frame(pi = c(0.1, 0.5, 0.4), mu = c(-3, 0, 3), sigma = sqrt(c(0.5, 
                                                                           0.75, 1)))
dat <- dat_generator(truth)

ggplot(data = dat %>% gather(key, value), aes(value)) + geom_density(aes(color = key)) + 
  theme_bw() + xlab("y") + ggtitle("y is mixture of {y1,y2,y3}")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

y <- dat$y
C <- 10  # to ensure large enough
N <- length(y)
input_dat <- list(y = y, N = N, C = C)
                                        # model_object<-stan_model(model_code=stan_model)
fit <- stan(file = "GaussianDP.stan", data = input_dat, iter = 1000, chains = 1)
results <- rstan::extract(fit)
