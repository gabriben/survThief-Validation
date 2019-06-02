load("../../meta.RData")

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
