## Comparison of rates of convergence (RMSE) of model averaging
## versus nonparametric kernel regression (kernel rate is -2/(4+p) where p is 
## the number of numeric predictors (-0.33 for p=2), correctly specified 
## parametric model's rate would be -0.50)

library(ma)
library(np)
options(np.messages=FALSE,np.tree=TRUE)
require(robustbase)

set.seed(42)
M <- 1000
N <- seq(500,1500,by=100)
nobs <- numeric()
mse.ma <- numeric()
mse.np <- numeric()
i <- 1
for(m in 1:M) {
    cat(paste("\rReplication ",m," of ",M,sep=""))
    for(n in N) {
        x1 <- runif(n,-1,1)
        x2 <- runif(n,-1,1)
        dgp <- sin(pi*x1) + x2**2 + x1*x2     
        y <- dgp + rnorm(n,sd=0.5*sd(dgp))
        mse.ma[i] <- mean((fitted(lm.ma(y~x1+x2,verbose=FALSE))-dgp)^2)
        mse.np[i] <- mean((fitted(npreg(y~x1+x2,ckertype="epanechnikov"))-dgp)^2)
        nobs[i] <- n
        i <- i + 1
    }
    ## Regress log(RMSE) on log(n), slope is empirical rate of convergence
    model.ma <- ltsReg(log(sqrt(mse.ma))~log(nobs))
    model.np <- ltsReg(log(sqrt(mse.np))~log(nobs))
    ylim=c(min(log(sqrt(c(mse.ma,mse.np)))),max(log(sqrt(c(mse.ma,mse.np)))))
    plot(log(nobs),log(sqrt(mse.ma)),
         ylim=ylim,
         xlab="log n",
         ylab="log RMSE",
         cex=.5,
         col=1,
         sub=paste("Replication ",m," of ",M,", MA rate = ",formatC(coef(model.ma)[2],format="f",digits=2),
                   ", NP rate = ",formatC(coef(model.np)[2],format="f",digits=2),sep=""))
    points(log(nobs),log(sqrt(mse.np)),col=2,cex=.5,pch=4)
    abline(model.ma)
    abline(model.np,col=2,lty=2)
    legend("topright",c("MA","NP"),col=1:2,lty=1:2,bty="n")

}
