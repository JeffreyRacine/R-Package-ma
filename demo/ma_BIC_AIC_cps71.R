## Compare model averaging, BIC, and AIC model selection for hold-out wage data

set.seed(42)

library(ma)
library(MASS)
data(cps71)
set.seed(42)

n <- nrow(cps71)
## There are 1.4 million possible hold-out samples of size 3
## choose(nrow(cps71),3)
## [1] 1414910
n.eval <- 3
n.train <- n - n.eval
M <- 1000

pse.ma <- numeric()
pse.lm.AIC <- numeric()
pse.lm.BIC <- numeric()

## Formula for BIC and AIC model selection - polynomial of degree 0 to 10,
## the same is used for model average approach

formula.upper <- formula(logwage~age + I(age^2) + I(age^3) + I(age^4) +
                         I(age^5) + I(age^6) + I(age^7) + I(age^8) +
                         I(age^9) + I(age^10))

for(m in 1:M) {

    cat(paste("\rReplication ",m," of ",M,sep=""))

    ## Split the sample into training and evaluation sets, compute
    ## the predicted square error on the independent evaluation data
    
    ii <- sample(seq(1,nrow(cps71)),replace=FALSE)

    cps71.train <- cps71[ii[1:n.train],]
    cps71.eval <- cps71[ii[(n.train+1):n],]

    ## Linear regression model, quadratic in age

    ## BIC model selection
    
    model.lm.AIC <- stepAIC(lm(logwage ~ age, data=cps71.train),
                            scope=list(upper=formula.upper,lower=~1),
                            k=2,
                            trace=0)
    
    pse.lm.AIC[m] <- mean((predict(model.lm.AIC,newdata=cps71.eval)-cps71.eval$logwage)^2)

    ## BIC model selection
    
    model.lm.BIC <- stepAIC(lm(logwage ~ age, data=cps71.train),
                              scope=list(upper=formula.upper,lower=~1),
                              k=log(length(cps71.train$logwage)),
                              trace=0)

    pse.lm.BIC[m] <- mean((predict(model.lm.BIC,newdata=cps71.eval)-cps71.eval$logwage)^2)

    model.ma <- lm.ma(logwage ~ age,
                      data= cps71.train,
                      degree.max = 10,
                      verbose=FALSE)

    pse.ma[m] <- mean((predict(model.ma,newdata=cps71.eval)-cps71.eval$logwage)^2)
    
    ## Create a boxplot of pse in real time and report the median relative efficiencies

    foo <- data.frame(pse.ma,pse.lm.BIC,pse.lm.AIC)
    names(foo) <- c("MA","BIC","AIC")
    foobar <- apply(cbind(pse.ma,pse.lm.BIC,pse.lm.AIC),2,median)
    foobar <- formatC(foobar/min(foobar),format="f",digits=3)
    boxplot(foo, outline=FALSE, notch=TRUE, ylab="PSE",
            sub=paste("Replication ",m," of ",M,", relative efficiency = ", 
                      foobar[1], ", ",foobar[2],", ",foobar[3],sep=""))    

}

