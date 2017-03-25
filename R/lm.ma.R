lm.ma <- function(...) UseMethod("lm.ma")

lm.ma.default <- function(y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("glp","tensor","additive"),
                          compute.deriv=FALSE,
                          p.max=NULL,
                          method=c("jma","mma"),
                          ma.weights=NULL,
                          bootstrap.ci=FALSE,
                          B=100,
                          alpha=0.05) {

    basis <- match.arg(basis)
    method <- match.arg(method)

    if(!is.null(p.max)) if(p.max < 2) stop("You must average over at least two models")
    if(is.null(y) | is.null(X)) stop("You must provide data for y and X")
    if(!is.null(X) & !is.null(X.eval) & NCOL(X)!=NCOL(X.eval)) stop("X and X.eval must contain the same number of predictors")

    if(is.null(p.max)) p.max <- round((10/NCOL(X))*(NROW(X)/100)^0.25)
    
    Est <- lm.ma.Est(y=y,
                     X=X,
                     X.eval=X.eval,
                     basis=basis,
                     compute.deriv=compute.deriv,
                     p.max=p.max,
                     method=method,
                     ma.weights=ma.weights)
    
    if(bootstrap.ci) {
    
        B <- 100
        if(is.null(X.eval)) {
            n <- NROW(X) 
        } else {
            n <- NROW(X.eval)
        }
    
        boot.mat <- matrix(NA,n,B)
        if(compute.deriv) boot.deriv.array <- array(NA,c(n,B,NCOL(X)))
    
        for(b in 1:B) {
            ii <- sample(1:n,replace=TRUE)
            out.boot <- lm.ma.Est(y=y[ii],
                                  X=X[ii,],
                                  X.eval=if(is.null(X.eval)) { X } else {X.eval},
                                  basis=basis,
                                  compute.deriv=compute.deriv,
                                  p.max=p.max,
                                  method=method,
                                  ma.weights=Est$ma.weights)
            boot.mat[,b] <- out.boot$fitted
            if(compute.deriv) for(k in 1:NCOL(X)) boot.deriv.array[,b,k] <- out.boot$deriv[,k]
        }

        Est$fitted.ci.l <- apply(boot.mat,1,quantile,prob=alpha/2,type=1)
        Est$fitted.ci.u <- apply(boot.mat,1,quantile,prob=1-alpha/2,type=1)
        
        if(compute.deriv) {
            Est$deriv.ci.l <- Est$deriv.ci.u <- matrix(NA,n,NCOL(X))
            for(k in 1:NCOL(X)) {
                Est$deriv.ci.l[,k] <- apply(boot.deriv.array,3,quantile,prob=alpha/2,type=1)
                Est$deriv.ci.u[,k] <- apply(boot.deriv.array,3,quantile,prob=1-alpha/2,type=1)        
            }
        }
    }

    Est$p.max <- p.max

    class(Est) <- "lm.ma"
    return(Est)

}

lm.ma.Est <- function(y=NULL,
                      X=NULL,
                      X.eval=NULL,
                      basis=c("glp","tensor","additive"),
                      compute.deriv=FALSE,
                      p.max=NULL,
                      method=c("jma","mma"),
                      ma.weights=NULL) {
    
    ## Divide into factors and numeric

    xztmp <- splitFrame(data.frame(X))
    x <- xztmp$x
    z <- xztmp$z
    if(!is.null(X.eval)) {
        xztmp <- splitFrame(data.frame(X.eval))
        xeval <- xztmp$x
        zeval <- xztmp$z
    } else {
        xeval <- NULL
        zeval <- NULL
    }
    rm(xztmp)

    deriv <- NULL
    if(compute.deriv) {
        deriv.mat <- array(NA,c(if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},p.max,NCOL(X)))
        deriv <- matrix(NA,if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},NCOL(X))
    }

    K <- numeric(length=p.max)
    ma.mat <- matrix(NA,NROW(X),p.max)
    fitted.mat <- matrix(NA,if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},p.max)

    for(p in 1:p.max) {
        
        if(is.null(ma.weights)) {
            
            model.ma <- suppressWarnings(crs:::predict.factor.spline(x,
                                                                     y,
                                                                     z,
                                                                     K=cbind(rep(p,NCOL(x)),rep(1,NCOL(x))),
                                                                     basis=basis)$model)
            
            K[p] <- model.ma$rank

            ma.mat[,p] <- Dmat.func(model.ma,method=method)
        
        }

        fitted.mat[,p] <- suppressWarnings(crs:::predict.factor.spline(x,
                                                                       y,
                                                                       z,
                                                                       xeval=xeval,
                                                                       zeval=zeval,
                                                                       K=cbind(rep(p,NCOL(x)),rep(1,NCOL(x))),
                                                                       basis=basis)$fitted.values[,1])
        
        if(compute.deriv) {
            for(k in 1:NCOL(X)) {
                model.deriv <- suppressWarnings(crs:::deriv.factor.spline(x,
                                                                          y,
                                                                          z,
                                                                          K=cbind(rep(p,NCOL(x)),rep(1,NCOL(x))),
                                                                          xeval=xeval,
                                                                          zeval=zeval,
                                                                          knots="uniform",
                                                                          basis=basis,
                                                                          deriv.index=k,
                                                                          deriv=1)[,1])
        
                deriv.mat[,p,k] <- model.deriv
            }
        }

    }

    if(is.null(ma.weights)) {

        ## Largest dimensional model is the last
        
        sigsq <- summary(model.ma)$sigma^2
        
        ## Solve the quadratic program for the Mallows model average weights
        
        M <- ncol(ma.mat)
        D <- t(ma.mat)%*%ma.mat
        if(qr(D)$rank<M) D <- D + diag(1e-10,M,M)
        A <- cbind(rep(1,M),diag(1,M,M))
        b0 <- c(1,rep(0,M))
        if(method=="mma") {
            d <- -sigsq*K
        } else {
            d <- t(y)%*%ma.mat
        }        
        b <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution

    } else {
        ## For bootstrapping use weights from initial call
        b <- ma.weights
    }

    if(compute.deriv) for(k in 1:NCOL(X)) deriv[,k] <- deriv.mat[,,k]%*%b

    return(list(fitted.values=fitted.mat%*%b,
                deriv=deriv,
                ma.weights=b,
                y=y,
                X=X,
                basis=basis,
                compute.deriv=compute.deriv,
                p.max=p.max,
                method=method))

    
}

lm.ma.formula <- function(formula,
                          data=list(),
                          y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("glp","tensor","additive"),
                          compute.deriv=FALSE,
                          p.max=NULL,
                          method=c("jma","mma"),
                          ma.weights=NULL,
                          bootstrap.ci=FALSE,
                          B=100,
                          alpha=0.05,
                          ...) {


  mf <- model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  tydat <- model.response(mf)
  txdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

  Est <- lm.ma.default(y=tydat,
                       X=txdat,
                       X.eval=NULL,
                       basis=basis,
                       compute.deriv=compute.deriv,
                       p.max=p.max,
                       method=method,
                       ma.weights=ma.weights,
                       bootstrap.ci=bootstrap.ci,
                       B=B,
                       alpha=alpha)

  Est$r.squared <- RSQfunc(tydat,Est$fitted.values)
  Est$residuals <- tydat - Est$fitted.values

  Est$call <- match.call()
  Est$formula <- formula
  Est$terms <- mt

  return(Est)

}

print.lm.ma <- function(x,
                        ...) {
    cat("Call:\n")
    print(x$call)

}

summary.lm.ma <- function(object,
                          ...) {

  cat("Call:\n")
  print(object$call)
  cat("\nModel Averaging Regression\n",sep="")

  cat("\n\n")

}

## Method for predicting given a new data frame.

predict.lm.ma <- function(object,
                          newdata=NULL,
                          ...) {

  if(is.null(newdata)) {
      fitted.values <- fitted(object)
  } else{
    tt <- terms(object)
    exdat <- model.frame(delete.response(tt),newdata,xlev=object$xlevels)
    Est <- lm.ma.default(y=object$y,
                         X=object$X,
                         X.eval=exdat,
                         basis=object$basis,
                         compute.deriv=object$compute.deriv,
                         p.max=object$p.max,
                         method=object$method,
                         ma.weights=object$ma.weights)
    fitted.values <- Est$fitted.values
    deriv <- Est$deriv
  }

  attr(fitted.values, "deriv") <- deriv

  return(fitted.values)

}
