library(survThief)

atRisk <- data.frame(time= c(0, 12 , 24, 36, 48, 60, 72, 84, 96, 108),
                     blue=c(100, 69, 60, 50, 31, 23, 16, 15, 10, 10),
                     orange=c(47, 33, 29, 27, 18, 18, 18, 18, 18, 18))

load("clickData.RData")

## clickData <- digi("Benrashid.jpg")

img <- survThief("Benrashid original.jpg", atRisk,
                 clickData = clickData, quantization = T)

## alternatively, redefine clicks:
img <- survThief("Benrashid original.jpg", atRisk, quantization = T)


## names(img)[1:2] <- c("hybrid", "open")
## img$survData$stratum[img$survData$stratum == "blue"] <- "hybrid"
## img$survData$stratum[img$survData$stratum == "orange"] <- "open"
save(img, file =  "Benrashid.RData")

load("Benrashid.RData")

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=img$survData)

library(survminer)
## pdf("Benrashid.pdf",width=14,height=10)
ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
           ## palette = substring(names(fit$strata), nchar("stratum= ")),
           risk.table.col = "strata", RUEpval = TRUE, 
           ggtheme = theme_bw(), risk.table.y.text=F,
           break.time.by=12)
## dev.off()
