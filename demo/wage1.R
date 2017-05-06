## Simple illustration with two numeric and two categorical predictors

library(ma)

data(wage1)
attach(wage1)
model <- lm.ma(lwage ~ female + married + educ + exper,
               compute.deriv=TRUE,
               parallel=TRUE)

apply(coef(model),2,summary)
plot(model)
plot(model,plot.deriv=TRUE)

summary(model)



