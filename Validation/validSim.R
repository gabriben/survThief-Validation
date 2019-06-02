
## @b : beta b
## @n : sample size per arm
## @a : number of arms
## @r : censoring rate

survSim <- function(b, n, a, r){

  sim <- array(numeric(), c(n,a))
    
  for(j in 1:a){

    m <- 100 # mixture amount
    B <- rbeta(m, 1, b)
    p <- numeric(m)
    p[1] <- B[1]
    p[2:m] <- sapply(2:m, function(i) B[i] * prod(1 - B[1:(i-1)]))

    alpha <- runif(m, 0.5, 1) # shape
    lambda <- runif(m, 5, 10) # scale

    y <- array(0,c(m,n))

    for(i in 1:n){
      y[i,] <- rweibull(n, alpha[i], lambda[i])
    }

    sim[,j] <- (p %*% y)[1,]
  }

  ## censoring
  C <- rexp(n * a, r)
  simC <- pmin(c(sim), C)
  
  sim <- cbind(simC, sort(rep(1:a, n)), as.numeric(c(sim) <= C))

  d <- data.frame(sim[order(sim[,1]),])
  colnames(d) <- c("time", "stratum", "status")

  return(d)

}

setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Validation")

library(survival)
library(survminer)
library(profvis)


## profvis({

sims <- 50
dims <- c(500,1000,2000)
dimsNames <- c("small", "medium", "large")
IPD <- rep(list(list()), 3)
atRisk <- rep(list(list()), 3)
names(atRisk) <- names(IPD) <- dimsNames

for(w in 1:3){
  for(s in 1:sims){
    d <- survSim(2, 25, 2, 0.1)

    ## IPD[[(w-1) * sims + s]] <- d
    IPD[[w]][[s]] <- d
    
    if(s<10) n <- paste0(dimsNames[w], "_sim", "0", s)
    else n <- paste0(dimsNames[w], "_sim", s)

    ## library(plotwidgets)
    library(colorspace)
    
    H <- sample(6,2) # sample without replacement, to sample from different parent colours
    H <- H * runif(2,0,2)
    H <- (360 * H/12)
    S <- runif(2,0.3,1)
    V <- runif(2,0.5,1)
    ## L <- runif(2,0.5,1)
    ## clrs <- hsl2rgb(t(cbind(H,S,L)))
    clrs <- (as(HSV(H,S,V), "RGB"))@coords

    ## plot(rep(1,2), 1:2, col = rgb(clrs), cex=10, axes = F, pch = 16)
    
    fit <- survfit(Surv(time, status) ~ stratum, data=d)
    p <- ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
                    palette = rgb(clrs),                      
                    risk.table.col = "strata", RUEpval = TRUE,
                    ggtheme = theme_bw(), risk.table.y.text=F)

    png(paste0("plots/", n, ".png"), width=dims[w], height=dims[w])
    print({p})
    dev.off()
    ## save(d, file = paste0(n, "/", n, ".RData"))


    temp <- p$data.survtable %>%
      .[names(.) %in% c("strata", "time", "n.risk")]

    T <- nrow(temp) / 2

    ## [[(w-1) * sims + s]]
    atRisk[[w]][[s]] <- data.frame(time = temp$time[1:T],
                                   treat1 = temp$n.risk[1:T],
                                   treat2 = temp$n.risk[(T+1):(2*T)])

    print(paste0(n, ".png"))
  }
}

save(atRisk, file = "atRisk.RData")
save(IPD, file = "IPD.RData")

## })
