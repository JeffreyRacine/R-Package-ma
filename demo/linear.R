## Demonstration of how one would replicate simple linear regression

require(ma)

set.seed(42)
n <- 1000

x1 <- rnorm(n)
x2 <- rnorm(n)
z1 <- rbinom(n,1,.5)
z2 <- rbinom(n,1,.5)

dgp <- x1 + x2 + z1 + z2

y <- dgp + rnorm(n,sd=0.5*sd(dgp))

z1 <- factor(z1)
z2 <- factor(z2)

model.lm <- lm(y~x1+x2+z1+z2)
model.ma <- lm.ma(y~x1+x2+z1+z2,
                  degree.min=1,
                  degree.max=1,
                  basis="additive",
                  compute.deriv=TRUE,
                  vc=FALSE)

summary(model.lm)
summary(model.ma)

coef(model.lm)[-1]
(coef(model.ma)[z1==1&z2==1,])[1,]
