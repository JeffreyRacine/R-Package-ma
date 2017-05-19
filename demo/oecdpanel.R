## Simple illustration with four numeric and two categorical predictors

library(ma)

data(oecdpanel)
attach(oecdpanel)

oecd <- factor(oecd)
year <- ordered(year)

model.ma <- lm.ma(growth ~ oecd + year + initgdp + popgro + inv + humancap)

summary(model.ma)

plot(model.ma,plot.data=TRUE,plot.rug=TRUE)
