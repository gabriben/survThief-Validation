simulWeib <- function(N, lambda, alpha, beta, rateC)
{
  ## covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))

  ## Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(x * beta)))^(1 / alpha)

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

d <- simulWeib(N=10, lambda=0.01, alpha=1, beta=-0.6, rateC=0.005)
fit <- survfit(Surv(time, status) ~ x, data=d)

png(paste0("sim.png"), width=1000, height=1000)
print({
  ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
             palette = substring(names(fit$strata), nchar("stratum= ")),
             risk.table.col = "strata", RUEpval = TRUE, 
             ggtheme = theme_bw(), risk.table.y.text=F)
})
dev.off()

save(d, file = paste0("sim.RData"))
