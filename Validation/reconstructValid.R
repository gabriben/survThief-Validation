detach("package:survThief", unload=TRUE)
install.packages("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Digitize Survival Curves in sections/Code/survThief", repos = NULL, type="source")
library(survThief)

atRisk <- data.frame(time= c(0, 12 , 24, 36, 48, 60, 72, 84, 96, 108),
                     blue=c(100, 69, 60, 50, 31, 23, 16, 15, 10, 10),
                     orange=c(47, 33, 29, 27, 18, 18, 18, 18, 18, 18))

for(i in 1:sims){
  d <- survSim(2, 25, 2, 0.1)
  if(i<10) n <- paste0("sim", "0", i)
  else n <- paste0("sim", i)

  

  fit <- survfit(Surv(time, status) ~ stratum, data=d)
  png(paste0(n, "/", n, ".png"), width=1000, height=1000)

  save(d, file = paste0(n, "/", n, ".RData"))
  i
}

