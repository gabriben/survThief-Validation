detach("package:survThief", unload=TRUE)
install.packages("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Digitize Survival Curves in sections/Code/survThief", repos = NULL, type="source")
library(survThief)

atRisk <- data.frame(time= c(0, 24 , 48, 72, 96, 120),
                     red=c(97, 55, 48, 35, 19, 12),
                     black=c(180, 102, 62, 33, 26, 14))


setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Reconstructed/Shrestha")

pic <- quantize("Shrestha no dots.png")
clickData <- digi(as.cimg(pic[,,, -4]))
## load("clickData.RData")

img <- survThief("Shrestha no dots.png", atRisk, quantization = T, method = "exact")

load("Shrestha.RData")
names(img)[1:2] <- c("open", "hybrid")
img$survData$stratum[img$survData$stratum == "black"] <- "hybrid"
img$survData$stratum[img$survData$stratum == "red"] <- "open"
save(img, file = "Shrestha.RData")

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=img$survData)

library(survminer)
## pdf("results/detectCensoring V2.pdf",width=14,height=10)
ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
           palette = substring(names(fit$strata), nchar("stratum= ")),
           risk.table.col = "strata", RUEpval = TRUE, 
           ggtheme = theme_bw(), risk.table.y.text=F,
           break.time.by=24)
## dev.off()

## save(demo, file="results/detectCensoring V3.RData")
