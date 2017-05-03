## Demonstration of anova-based significance test

library(ma)

n <- 500
x1 <- runif(n)
x2 <- runif(n)

## One irrelevant predictor (x2)

dgp <- sin(2*pi*x1)
y <- dgp + rnorm(n,sd=0.25*sd(dgp))
model <- lm.ma(y~x1+x2,compute.anova=TRUE)
summary(model)

## Both predictors relevant

dgp <- sin(2*pi*x1)+cos(2*pi*x2)
y <- dgp + rnorm(n,sd=0.25*sd(dgp))
model <- lm.ma(y~x1+x2,compute.anova=TRUE)
summary(model)

## Both predictors irrelevant

y <- rnorm(n,sd=0.25*sd(dgp))
model <- lm.ma(y~x1+x2,compute.anova=TRUE)
summary(model)

