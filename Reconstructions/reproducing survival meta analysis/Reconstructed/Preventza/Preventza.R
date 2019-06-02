detach("package:survThief", unload=TRUE)
install.packages("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Digitize Survival Curves in sections/Code/survThief", repos = NULL, type="source")
library(survThief)

atRisk <- data.frame(time= c(0, 1 , 2, 3, 4, 5),
                     blue=c(45, 27, 24, 18, 8, 3),
                     red=c(274, 194, 163, 140, 112, 88))


setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Reconstructed/Preventza")

## load("clickData.RData")

img <- survThief("Preventza original.png", atRisk, quantization = T, method = "reconstruct")

names(img)[1:2] <- c("hybrid", "open")
img$survData$stratum[img$survData$stratum == "blue"] <- "hybrid"
img$survData$stratum[img$survData$stratum == "red"] <- "open"
save(img, file =  "Preventza.RData")

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=img$survData)

library(survminer)
## pdf("results/detectCensoring V2.pdf",width=14,height=10)
ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
           ## palette = substring(names(fit$strata), nchar("stratum= ")),
           risk.table.col = "strata", RUEpval = TRUE, 
           ggtheme = theme_bw(), risk.table.y.text=F,
           break.time.by=1)
## dev.off()

## save(demo, file="results/detectCensoring V3.RData")