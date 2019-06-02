setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Reconstructions/reproducing survival meta analysis/Reconstructed")

studies <- list.dirs(getwd(), recursive = F, full.names = F)

## each dataset is saved in an object "img" in .RData format.
## loading RData does not allow for reassigning objects to other names, so img is kept all along.
load(paste0(studies[1], "/", studies[1], ".RData"))
d <- cbind(img$survData, study = studies[1])
for(i in studies[-1]){
  load(paste0(i, "/", i, ".RData"))
  d <- rbind(d, cbind(img$survData, study = i))
}
d <- d[order(d$time),]

d$medAge <- 65 * (d$study == "Benrashid") * (d$stratum == "hybrid") +
  55 * (d$study == "Benrashid") * (d$stratum == "open") +
  63 * (d$study == "Shrestha") * (d$stratum == "hybrid") +
  62 * (d$study == "Shrestha") * (d$stratum == "open") +
  63 * (d$study == "Preventza") * (d$stratum == "hybrid") +
  68 * (d$study == "Preventza") * (d$stratum == "open") +
  74.7 * (d$study == "Yoshitake") * (d$stratum == "hybrid") +
  67.8 * (d$study == "Yoshitake") * (d$stratum == "open")

setwd("../")
save(d, file = "meta.RData")
