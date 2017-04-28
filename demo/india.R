## Estimate a semiparametric additive model averaged model

library(ma)

data(india)
attach(india)

faccsex <- factor(csex)
facctwin <- factor(ctwin)
faccbirthorder <- factor(cbirthorder)
facmunemployed <- factor(munemployed)
facmreligion <- factor(mreligion)
faccar <- factor(car)

model <- lm.ma(cheight ~ faccsex + facctwin + faccbirthorder + facmunemployed
               + facmreligion + faccar + cage + mbmi + medu,
               basis="additive",
               vc=FALSE)

summary(model)

plot(model,plot.data=TRUE)
plot(model,plot.deriv=TRUE)
