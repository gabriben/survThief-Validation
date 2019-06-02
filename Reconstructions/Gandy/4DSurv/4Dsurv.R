detach("package:survThief", unload=TRUE)
install.packages("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Digitize Survival Curves in sections/Code/survThief", repos = NULL, type="source")
library(survThief)

atRisk <- data.frame(time= c(0, 3 , 6, 9, 12),
                     blue=c(151, 47, 22, 10, 0),
                     red=c(151, 53, 28, 12, 0))


setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/Gandy/4DSurv")

## load("clickData.RData")

## clickData <- digi("Benrashid.jpg")

img <- survThief("4DSurv.png", atRisk, quantization = T)

names(img)[1:2] <- c("low risk", "high risk")
img$survData$stratum[img$survData$stratum == "blue"] <- "low risk"
img$survData$stratum[img$survData$stratum == "red"] <- "high risk"
save(img, file =  "4DSurv.RData")

load("Benrashid.RData")

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=img$survData)

library(survminer)
## pdf("testGapStraight.pdf",width=14,height=10)
ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
           palette = substring(names(fit$strata), nchar("stratum= ")),
           risk.table.col = "strata", RUEpval = TRUE, 
           ggtheme = theme_bw(), risk.table.y.text=F,
           break.time.by=3)
## dev.off()

## save(demo, file="results/detectCensoring V3.RData")
