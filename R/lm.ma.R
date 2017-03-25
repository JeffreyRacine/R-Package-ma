lm.ma <- function(...) UseMethod("lm.ma")

lm.ma.default <- function(y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("glp","tensor","additive"),
                          compute.deriv=TRUE,
                          p.max=25,
                          method=c("jma","mma"),
                          ma.weights=NULL) {

    basis <- match.arg(basis)
    method <- match.arg(method)

    if(p.max < 2) stop("You must average over at least two models")
    if(is.null(y) | is.null(X)) stop("You must provide data for y and X")
    if(!is.null(X) & !is.null(X.eval) & NCOL(X)!=NCOL(X.eval)) stop("X and X.eval must contain the same number of predictors")

    Est <- lm.ma.Est(y=y,
                     X=X,
                     X.eval=X.eval,
                     basis=basis,
                     compute.deriv=compute.deriv,
                     p.max=p.max,
                     method=method,
                     ma.weights=ma.weights)

    class(Est) <- "lm.ma"
    return(Est)

}

lm.ma.Est <- function(y=NULL,
                      X=NULL,
                      X.eval=NULL,
                      basis=c("glp","tensor","additive"),
                      compute.deriv=TRUE,
                      p.max=25,
                      method=c("jma","mma"),
                      ma.weights=NULL) {
    
    ## Divide into factors and numeric, if one regressor must be numeric

#    if(NCOL(X) > 1) {
        xztmp <- crs:::splitFrame(data.frame(X))
        x <- xztmp$x
        z <- xztmp$z
        if(!is.null(X.eval)) {
            xztmp <- crs:::splitFrame(data.frame(X.eval))
            xeval <- xztmp$x
            zeval <- xztmp$z
        } else {
            xeval <- NULL
            zeval <- NULL
        }
        rm(xztmp)
#    } else {
#        #if(!is.numeric(X)) stop("Single predictor must be of type numeric")
#        print(class(X))
#        stop()
#        x <- X
#        xeval <- X.eval
#        z <- NULL
#        zeval <- NULL
#    }

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
                          compute.deriv=TRUE,
                          p.max=25,
                          method=c("jma","mma"),
                          ma.weights=NULL,
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
                       ma.weights=ma.weights)

  Est$r.squared <- crs:::RSQfunc(tydat,Est$fitted.values)
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
#    has.ey <- succeedWithResponse(tt, newdata)
#    if (has.ey) {
#      eydat <- model.response(model.frame(tt,newdata))
#    } else {
#      eydat <- NULL
#    }
    exdat <- model.frame(delete.response(tt),newdata,xlev=object$xlevels)

    ## Return the predicted values.

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
