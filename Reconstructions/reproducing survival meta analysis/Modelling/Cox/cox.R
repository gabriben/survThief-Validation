
load("../../meta.RData")

fits$cox <- coxph(Surv(time, status) ~ stratum + study + medAge, data = d, ties="efron")
names(d)

ggforest(fits$cox)
