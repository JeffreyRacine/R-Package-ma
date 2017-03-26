lm.ma <- function(...) UseMethod("lm.ma")

lm.ma.default <- function(y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("glp","tensor","additive"),
                          compute.deriv=FALSE,
                          deriv.order=1,
                          degree.max=NULL,
                          knots=FALSE,
                          S=10,
                          exhaustive=TRUE,
                          method=c("jma","mma"),
                          ma.weights=NULL,
                          bootstrap.ci=FALSE,
                          B=199,
                          alpha=0.05,
                          weights=NULL,
                          vc=TRUE,
                          verbose=TRUE,
                          ...) {

    basis <- match.arg(basis)
    method <- match.arg(method)
    
    if(verbose) {
        options(crs.messages=FALSE)
    } else {
        options(crs.messages=TRUE)
    }
    
    if(!is.null(degree.max)) if(degree.max < 2) stop("You must average over at least two models")
    if(is.null(y) | is.null(X)) stop("You must provide data for y and X")
    if(!is.null(X) & !is.null(X.eval) & NCOL(X)!=NCOL(X.eval)) stop("X and X.eval must contain the same number of predictors")

    Est <- lm.ma.Est(y=y,
                     X=X,
                     X.eval=X.eval,
                     basis=basis,
                     compute.deriv=compute.deriv,
                     deriv.order=deriv.order,
                     degree.max=degree.max,
                     knots=knots,
                     S=S,
                     exhaustive=exhaustive,
                     method=method,
                     ma.weights=ma.weights,
                     weights=weights,
                     vc=vc,
                     verbose=verbose)
    
    if(bootstrap.ci) {
    
        B <- 100
        if(is.null(X.eval)) {
            n <- NROW(X) 
        } else {
            n <- NROW(X.eval)
        }
    
        boot.mat <- matrix(NA,n,B)
        if(compute.deriv) boot.deriv.array <- array(NA,c(n,B,Est$num.x))
    
        for(b in 1:B) {
            ii <- sample(1:n,replace=TRUE)
            out.boot <- lm.ma.Est(y=y[ii],
                                  X=X[ii,],
                                  X.eval=if(is.null(X.eval)) { X } else {X.eval},
                                  basis=basis,
                                  compute.deriv=compute.deriv,
                                  deriv.order=deriv.order,
                                  degree.max=degree.max,
                                  knots=knots,
                                  S=S,
                                  exhaustive=exhaustive,
                                  method=method,
                                  ma.weights=Est$ma.weights,
                                  weights=weights,
                                  vc=vc,
                                  verbose=verbose)
            boot.mat[,b] <- out.boot$fitted
            if(compute.deriv) for(k in 1:Est$num.x) boot.deriv.array[,b,k] <- out.boot$deriv[,k]
        }

        Est$fitted.ci.l <- apply(boot.mat,1,quantile,prob=alpha/2,type=1)
        Est$fitted.ci.u <- apply(boot.mat,1,quantile,prob=1-alpha/2,type=1)
        
        if(compute.deriv) {
            Est$deriv.ci.l <- Est$deriv.ci.u <- matrix(NA,n,Est$num.x)
            for(k in 1:Est$num.x) {
                Est$deriv.ci.l[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=alpha/2,type=1)
                Est$deriv.ci.u[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=1-alpha/2,type=1)        
            }
        }
    }

    class(Est) <- "lm.ma"
    return(Est)

}

lm.ma.Est <- function(y=NULL,
                      X=NULL,
                      X.eval=NULL,
                      basis=c("glp","tensor","additive"),
                      compute.deriv=FALSE,
                      deriv.order=1,
                      degree.max=NULL,
                      knots=FALSE,
                      S=10,
                      exhaustive=TRUE,
                      method=c("jma","mma"),
                      ma.weights=NULL,
                      weights=NULL,
                      vc=TRUE,
                      verbose=TRUE) {
    
    ## Divide into factors and numeric
    if(!vc) {
        xztmp <- splitFrame(as.data.frame(X))
    } else {
        xztmp <- splitFrame(as.data.frame(X),factor.to.numeric=TRUE)
    }
    xnames <- xztmp$xnames
    znames <- xztmp$znames
    x <- xztmp$x
    z <- xztmp$z
    num.x <- xztmp$num.x
    num.z <- xztmp$num.z
    is.ordered.z <- xztmp$is.ordered.z
    if(!is.null(X.eval)) {
        if(!vc) {
            xztmp <- splitFrame(as.data.frame(X.eval))
        } else {
            xztmp <- splitFrame(as.data.frame(X.eval),factor.to.numeric=TRUE)
        }
        xeval <- xztmp$x
        zeval <- xztmp$z
    } else {
        xeval <- NULL
        zeval <- NULL
    }
    rm(xztmp)
    if(is.null(z)) {
        include <- NULL
        lambda <- NULL
    } else {
        include <- rep(1,num.z)
        lambda <- rep(sqrt(.Machine$double.eps),num.z)
    }

    if(is.null(degree.max)) {
        if(knots) S <- S/2
        degree.max <- max(2,round((S/num.x)*(NROW(X)/100)^0.25))
    }
    
    P <- degree.max
    
    if(exhaustive) {
        
        degree.min = deriv.order
        
        ## Exhaustive evaluation over all combinations of K, search over
        ## lambda for each combination
        if(knots) {
            K.mat <- crs:::matrix.combn(K.vec1=degree.min:degree.max, K.vec2=degree.min:degree.max,num.x=num.x)
        } else {
            K.mat <- crs:::matrix.combn(K.vec1=degree.min:degree.max,num.x=num.x)
            K.mat <- cbind(K.mat[,1:num.x],matrix(1,nrow(K.mat),num.x,byrow=TRUE))
        }

        P <- NROW(K.mat)

    } else {
        ## Non-exhaustive
        if(knots) {
            K.mat <- cbind(rep(degree.min:degree.max,length=(degree.max-degree.min+1)*num.x),rep(1,length=(degree.max-degree.min+1)*num.x))
        } else {
            K.mat <- cbind(rep(degree.min:degree.max,length=(degree.max-degree.min+1)*num.x),rep(1,length=(degree.max-degree.min+1)*num.x))
        }
    }
    
    deriv <- NULL
    if(compute.deriv) {
        deriv.mat <- array(NA,c(if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},P,num.x))
        deriv <- matrix(NA,if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},num.x)
        colnames(deriv) <- xnames
    }

    K.rank <- numeric(length=degree.max)
    ma.mat <- matrix(NA,NROW(X),P)
    fitted.mat <- matrix(NA,if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},P)

    for(p in 1:P) {
        
        if(!exhaustive) {
            if(knots) {
                DS <- cbind(rep(p+deriv.order-1,num.x),rep(p+deriv.order-1,num.x))
            } else {
                DS <- cbind(rep(p+deriv.order-1,num.x),rep(1,num.x))
            }
        } else {
            DS <- cbind(K.mat[p,1:num.x],K.mat[p,(num.x+1):(2*num.x)])   
        }
        
        if(verbose) cat(paste("\rCandidate model ",p," of ",P," (degree.max = ",degree.max,")...",sep=""))

        if(is.null(ma.weights)) {

            if(vc & !is.null(num.z)) {

                model.ma <- suppressWarnings(crs:::predict.kernel.spline(x,
                                                                         y,
                                                                         z,
                                                                         K=DS,
                                                                         lambda=lambda,
                                                                         is.ordered.z=is.ordered.z,
                                                                         basis=basis,
                                                                         weights=weights,
                                                                         model.return=TRUE)$model[[1]])
            } else {
            
                model.ma <- suppressWarnings(crs:::predict.factor.spline(x,
                                                                         y,
                                                                         z,
                                                                         K=DS,
                                                                         I=include,
                                                                         basis=basis,
                                                                         weights=weights)$model)

            }

            K.rank[p] <- model.ma$rank

            ma.mat[,p] <- Dmat.func(model.ma,method=method)
        
        }

        if(vc & !is.null(num.z)) {

            fitted.mat[,p] <- suppressWarnings(crs:::predict.kernel.spline(x,
                                                                           y,
                                                                           z,
                                                                           xeval=xeval,
                                                                           zeval=zeval,
                                                                           K=DS,
                                                                           lambda=lambda,
                                                                           is.ordered.z=is.ordered.z,
                                                                           basis=basis,
                                                                           weights=weights)$fitted.values[,1])
            
        } else {

            fitted.mat[,p] <- suppressWarnings(crs:::predict.factor.spline(x,
                                                                           y,
                                                                           z,
                                                                           xeval=xeval,
                                                                           zeval=zeval,
                                                                           K=DS,
                                                                           I=include,
                                                                           basis=basis,
                                                                           weights=weights)$fitted.values[,1])
            
        }
        
        if(compute.deriv) {
            for(k in 1:num.x) {

                if(vc & !is.null(num.z)) {

                    model.deriv <- suppressWarnings(crs:::deriv.kernel.spline(x,
                                                                              y,
                                                                              z,
                                                                              K=DS,
                                                                              lambda=lambda,
                                                                              is.ordered.z=is.ordered.z,
                                                                              xeval=xeval,
                                                                              zeval=zeval,
                                                                              basis=basis,
                                                                              deriv.index=k,
                                                                              deriv=deriv.order,
                                                                              weights=weights)[,1])
                } else {

                    model.deriv <- suppressWarnings(crs:::deriv.factor.spline(x,
                                                                              y,
                                                                              z,
                                                                              K=DS,
                                                                              I=include,
                                                                              xeval=xeval,
                                                                              zeval=zeval,
                                                                              basis=basis,
                                                                              deriv.index=k,
                                                                              deriv=deriv.order,
                                                                              weights=weights)[,1])

                }
        
                deriv.mat[,p,k] <- model.deriv
            }
        }

    }
    
    if(verbose) cat("\r                                                    ")
    if(verbose) cat("\rComputing model average weights...")

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
            d <- -sigsq*K.rank
        } else {
            d <- t(y)%*%ma.mat
        }        
        b <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution

    } else {
        ## For bootstrapping use weights from initial call
        b <- ma.weights
    }

    if(compute.deriv) {
        if(verbose) cat("\r                                                    ")
        if(verbose) cat("\rComputing derivatives...")
        for(k in 1:num.x) deriv[,k] <- deriv.mat[,,k]%*%b
    }

    if(verbose) cat("\r                                                    ")

    return(list(fitted.values=fitted.mat%*%b,
                deriv=deriv,
                ma.weights=abs(b),
                y=y,
                X=X,
                basis=basis,
                compute.deriv=compute.deriv,
                deriv.order=deriv.order,
                degree.max=degree.max,
                knots=knots,
                S=S,
                exhaustive=exhaustive,
                method=method,
                num.x=num.x,
                num.z=num.z,
                rank.vec=K.rank,
                nobs=NROW(y),
                DS=K.mat))

}

lm.ma.formula <- function(formula,
                          data=list(),
                          y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("glp","tensor","additive"),
                          compute.deriv=FALSE,
                          deriv.order=1,
                          degree.max=NULL,
                          knots=FALSE,
                          S=10,
                          exhaustive=TRUE,
                          method=c("jma","mma"),
                          ma.weights=NULL,
                          bootstrap.ci=FALSE,
                          B=199,
                          alpha=0.05,
                          weights=NULL,
                          vc=TRUE,
                          verbose=TRUE,
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
                       deriv.order=deriv.order,
                       degree.max=degree.max,
                       knots=knots,
                       S=S,
                       exhaustive=exhaustive,
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

coef.lm.ma <- function(object,
                       ...) {

    if(!is.null(object$deriv)) {
        object$deriv
    } else {
        cat("\nlm.ma(...) was called with option compute.deriv=FALSE\n")
    }

    }

summary.lm.ma <- function(object,
                          ...) {

  cat("Call:\n")
  print(object$call)
  cat("\nModel Averaging Linear Regression\n",sep="")
  cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4), sep=""))
  cat(paste("\nMaximum basis degree: ", object$degree.max, sep=""))  
  cat(paste("\nNumber of observations: ", object$nobs, sep=""))
  cat(paste("\nRank of model averaged model frame: ", round(sum(object$rank.vec*object$ma.weights)), sep=""))
  cat(paste("\nResidual standard error: ", format(sqrt(sum(object$residuals^2)/(object$nobs-sum(object$rank.vec*object$ma.weights))),digits=4),
                                                  " on ", format(round(object$nobs-sum(object$rank.vec*object$ma.weights)))," degrees of freedom",sep=""))
  cat(paste("\nModel average criterion: ", object$method, sep=""))
  cat("\nNon-zero model ranks: ")
  cat(object$rank.vec[object$ma.weights>1e-03])
  cat("\nNon-zero model average weights: ")
  cat(formatC(object$ma.weights[object$ma.weights>1e-03],format="f",digits=3))
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
                         deriv.order=object$deriv.order,
                         degree.max=object$degree.max,
                         knots=object$knots,
                         S=object$S,
                         exhaustive=object$exhaustive,
                         method=object$method,
                         ma.weights=object$ma.weights)
    fitted.values <- Est$fitted.values
    deriv <- Est$deriv
  }

  attr(fitted.values, "deriv") <- deriv

  if(is.null(Est$deriv)) {
      return(fitted.values)
  } else {
      return(list(fit=fitted.values,deriv=deriv))
  }

}
