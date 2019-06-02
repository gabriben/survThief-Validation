detach("package:survThief", unload=TRUE)
install.packages("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Digitize Survival Curves in sections/Code/survThief", repos = NULL, type="source")
library(survThief)

atRisk <- data.frame(time= c(0, 12 , 24, 36, 48, 60, 72, 84, 96, 108),
                     blue=c(100, 69, 60, 50, 31, 23, 16, 15, 10, 10),
                     orange=c(47, 33, 29, 27, 18, 18, 18, 18, 18, 18))


setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Reconstructed/Benrashid")

load("clickData.RData")

clickData <- digi("Benrashid.jpg")

img <- survThief("Benrashid original.jpg", atRisk, clickData = clickData,
                 quantization = T, stairs = F)

names(img)[1:2] <- c("hybrid", "open")
img$survData$stratum[img$survData$stratum == "blue"] <- "hybrid"
img$survData$stratum[img$survData$stratum == "orange"] <- "open"
save(img, file =  "Benrashid.RData")

load("Benrashid.RData")

library(survival)
fit <- survfit(Surv(time, status) ~ stratum, data=img$survData)

library(survminer)
## pdf("testGapStraight.pdf",width=14,height=10)
ggsurvplot(fit, risk.table = T, tables.theme = theme_cleantable(),
           ## palette = substring(names(fit$strata), nchar("stratum= ")),
           risk.table.col = "strata", RUEpval = TRUE, 
           ggtheme = theme_bw(), risk.table.y.text=F,
           break.time.by=12)
## dev.off()

## save(demo, file="results/detectCensoring V3.RData")

install.packages("image.LineSegmentDetector", repos = "https://bnosac.github.io/drat")
install.packages("pixmap")

library(image.LineSegmentDetector)
library(pixmap)
image <- read.pnm(file = "Benrashid original.jpg", cellres = 1)
plot(image)
linesegments <- image_line_segment_detector(image@grey * 255)
plot(linesegments)


library(magick)
f <- tempfile(fileext = ".pgm")
x <- image_read("Benrashid original.jpg")
x <- image_convert(x, format = "pgm", depth = 8)
image_write(x, path = f, format = "pgm")

image <- read.pnm(file = f, cellres = 1)
linesegments <- image_line_segment_detector(image@grey * 255)
plot(image)
plot(linesegments, add = TRUE, col = "red")
