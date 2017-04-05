## Function call.

lm.ma <- function(...) UseMethod("lm.ma")

## Default formula interface.

lm.ma.formula <- function(formula,
                          data=list(),
                          y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          alpha=0.05,
                          B=199,
                          basis.vec=NULL,
                          basis=c("auto","tensor","taylor","additive"),
                          bootstrap.ci=FALSE,
                          compute.anova=FALSE,
                          compute.deriv=FALSE,
                          degree.max=NULL,
                          degree.min=1,
                          deriv.index=NULL,
                          deriv.order=1,
                          K.mat=NULL,
                          knots=FALSE,
                          lambda=1e-02,
                          ma.weights=NULL,
                          method=c("jma","mma"),
                          rank.vec=NULL,
                          S=2,
                          segments.max=3,
                          tol=1e-12,
                          vc=TRUE,
                          verbose=TRUE,
                          weights=NULL,
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
                       degree.min=degree.min,
                       deriv.index=deriv.index,
                       degree.max=degree.max,
                       lambda=lambda,
                       segments.max=segments.max,
                       knots=knots,
                       S=S,
                       method=method,
                       ma.weights=ma.weights,
                       basis.vec=basis.vec,
                       rank.vec=rank.vec,
                       bootstrap.ci=bootstrap.ci,
                       B=B,
                       alpha=alpha,
                       weights=weights,
                       vc=vc,
                       verbose=verbose,
                       tol=tol,
                       compute.anova=compute.anova,
                       ...)
  
  Est$r.squared <- RSQfunc(tydat,Est$fitted.values)
  Est$residuals <- tydat - Est$fitted.values

  Est$call <- match.call()
  Est$formula <- formula
  Est$terms <- mt

  return(Est)

}

## Default function.

lm.ma.default <- function(y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          alpha=0.05,
                          B=199,
                          basis.vec=NULL,
                          basis=c("auto","tensor","taylor","additive"),
                          bootstrap.ci=FALSE,
                          compute.anova=FALSE,
                          compute.deriv=FALSE,
                          degree.max=NULL,
                          degree.min=1,
                          deriv.index=NULL,
                          deriv.order=1,
                          K.mat=NULL,
                          knots=FALSE,
                          lambda=1e-02,
                          ma.weights=NULL,
                          method=c("jma","mma"),
                          rank.vec=NULL,
                          S=2,
                          segments.max=3,
                          tol=1e-12,
                          vc=TRUE,
                          verbose=TRUE,
                          weights=NULL,
                          ...) {

    basis <- match.arg(basis)
    method <- match.arg(method)

    if(!is.null(degree.max)) if(degree.max < 2) stop("You must average over at least two models")
    if(is.null(y) | is.null(X)) stop("You must provide data for y and X")
    if(!is.null(X) & !is.null(X.eval) & NCOL(X)!=NCOL(X.eval)) stop("X and X.eval must contain the same number of predictors")
    if(compute.deriv & (degree.min < deriv.order)) stop("Minimum degree (degree.min) must be at least as large the order of the derivative required (deriv.order)")
    if(degree.min < 1) stop("Minimum degree (degree.min) must be at least one")
    if(!is.null(deriv.index)) if(any(deriv.index < 1) | any(deriv.index > NCOL(X))) stop("Derivative indices must correspond to columns of X")

    ## First obtain weights, then in subsequent call computes fits and
    ## derivatives

    Est <- lm.ma.Est(y=y,
                     X=X,
                     X.eval=X.eval,
                     basis=basis,
                     compute.deriv=FALSE,
                     deriv.order=deriv.order,
                     degree.min=degree.min,
                     deriv.index=deriv.index,
                     degree.max=degree.max,
                     lambda=lambda,
                     segments.max=segments.max,
                     knots=knots,
                     S=S,
                     method=method,
                     ma.weights=ma.weights,
                     basis.vec=basis.vec,
                     rank.vec=rank.vec,
                     K.mat=K.mat,
                     weights=weights,
                     vc=vc,
                     verbose=verbose,
                     tol=tol,
                     ...)
    
    ## Save rank vector and matrix of degrees and segments (not
    ## computed in next call since ma.weights is passed, but needed
    ## for summary)

    if(compute.deriv | !is.null(X.eval)) {
        rank.vec <- Est$rank.vec
        DS <- Est$DS
        basis.vec <- Est$basis.vec

        Est <- lm.ma.Est(y=y,
                         X=X,
                         X.eval=X.eval,
                         basis=basis,
                         compute.deriv=compute.deriv,
                         deriv.order=deriv.order,
                         degree.min=degree.min,
                         deriv.index=deriv.index,
                         degree.max=degree.max,
                         lambda=lambda,
                         segments.max=segments.max,
                         knots=knots,
                         S=S,
                         method=method,
                         ma.weights=Est$ma.weights,
                         basis.vec=Est$basis.vec,
                         rank.vec=Est$rank.vec,
                         K.mat=Est$DS,
                         weights=weights,
                         vc=vc,
                         verbose=verbose,
                         tol=tol,
                         ...)
        
        Est$rank.vec <- rank.vec
        Est$DS <- DS
        Est$basis.vec <- basis.vec

    }

    ## Bootstrap if requested uses pre-computed ma.weights
    
    if(bootstrap.ci) {

        if(is.null(X.eval)) {
            n <- NROW(X) 
        } else {
            n <- NROW(X.eval)
        }
        
        ncol.X <- NCOL(X)
    
        boot.mat <- matrix(NA,n,B)

        if(is.null(deriv.index)) deriv.index <- 1:ncol.X
        num.deriv <- length(deriv.index)

        if(compute.deriv) boot.deriv.array <- array(NA,c(n,B,num.deriv))

        for(b in 1:B) {
            if(verbose) cat(paste("\rBootstrap replication ",b," of ",B,sep=""))
            ii <- sample(1:n,replace=TRUE)
            out.boot <- lm.ma.Est(y=y[ii],
                                  X=X[ii,],
                                  X.eval=if(is.null(X.eval)){X}else{X.eval},
                                  basis=basis,
                                  compute.deriv=compute.deriv,
                                  deriv.order=deriv.order,
                                  degree.min=degree.min,
                                  deriv.index=deriv.index,
                                  degree.max=degree.max,
                                  lambda=lambda,
                                  segments.max=segments.max,
                                  knots=knots,
                                  S=S,
                                  method=method,
                                  ma.weights=Est$ma.weights,
                                  basis.vec=Est$basis.vec,
                                  rank.vec=Est$rank.vec,
                                  K.mat=Est$DS,
                                  weights=weights,
                                  vc=vc,
                                  verbose=FALSE,
                                  tol=tol,
                                  ...)
            boot.mat[,b] <- out.boot$fitted
            if(compute.deriv) for(k in 1:num.deriv) boot.deriv.array[,b,k] <- out.boot$deriv[,k]
        }

        if(verbose) cat("\r                                     ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing quantiles...")

        Est$fitted.ci.l <- apply(boot.mat,1,quantile,prob=alpha/2,type=1,na.rm=TRUE)
        Est$fitted.ci.u <- apply(boot.mat,1,quantile,prob=1-alpha/2,type=1,na.rm=TRUE)
        Est$fitted.scale <- apply(boot.mat,1,mad,na.rm=TRUE)
        
        if(compute.deriv) {
            Est$deriv.ci.l <- Est$deriv.ci.u <- Est$deriv.scale <- matrix(NA,n,num.deriv)
            for(k in 1:num.deriv) {
                Est$deriv.ci.l[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=alpha/2,type=1,na.rm=TRUE)
                Est$deriv.ci.u[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=1-alpha/2,type=1,na.rm=TRUE)     
                Est$deriv.scale[,k] <- apply(boot.deriv.array[,,k],1,mad,na.rm=TRUE)  
            }
        }
        if(verbose) cat("\r                                     ")
        if(verbose) cat("\r")
    }

    if(compute.anova) {

        nrow.X <- NROW(X)
        ncol.X <- NCOL(X)
        
        P.vec <- numeric(length=ncol.X)
        F.stat <- numeric(length=ncol.X)
        F.boot <- numeric(length=B)
        
        if(!is.null(X.eval)) {
            Est.ssu <- lm.ma.Est(y=y,
                                 X=X,
                                 X.eval=NULL,
                                 basis=basis,
                                 compute.deriv=FALSE,
                                 deriv.order=deriv.order,
                                 degree.min=degree.min,
                                 deriv.index=deriv.index,
                                 degree.max=degree.max,
                                 lambda=lambda,
                                 segments.max=segments.max,
                                 knots=knots,
                                 S=S,
                                 method=method,
                                 ma.weights=Est$ma.weights,
                                 basis.vec=Est$basis.vec,
                                 rank.vec=Est$rank.vec,
                                 K.mat=Est$DS,
                                 weights=weights,
                                 vc=vc,
                                 verbose=FALSE,
                                 tol=tol,
                                 ...)
        } else {
            Est.ssu <- Est
        }

        ssu <- sum((y-Est.ssu$fitted.values)^2) 
        ssu.rank <- Est.ssu$ma.model.rank

        for(k in 1:ncol.X) {

            if(verbose) cat(paste("\rAnova for predictor ",k," of ",ncol.X,sep=""))

            is.numeric.X.k <- is.numeric(X[,k])
            
            ## With > 1 predictor, restricted model does not incorporate the kth predictor

            if(ncol.X>1 & (Est.ssu$num.x>1 | (Est.ssu$num.x==1 & !is.numeric.X.k))) {
                X.res <- X[,-k,drop=FALSE]
                Est.ssr <- lm.ma.Est(y=y,
                                     X=X.res,
                                     X.eval=NULL,
                                     basis=basis,
                                     compute.deriv=FALSE,
                                     deriv.order=deriv.order,
                                     degree.min=degree.min,
                                     deriv.index=deriv.index,
                                     degree.max=degree.max,
                                     lambda=lambda,
                                     segments.max=segments.max,
                                     knots=knots,
                                     S=S,
                                     method=method,
                                     ma.weights=Est.ssu$ma.weights,
                                     basis.vec=Est.ssu$basis.vec,
                                     rank.vec=Est.ssu$rank.vec,
                                     K.mat=if(is.numeric.X.k){Est.ssu$DS[,c(-k,-(k+Est.ssu$num.x)),drop=FALSE]}else{Est.ssu$DS},
                                     weights=weights,
                                     vc=vc,
                                     verbose=FALSE,
                                     tol=tol,
                                     ...)

                ssr <- sum((y-Est.ssr$fitted.values)^2) 
                ssr.rank <- Est.ssr$ma.model.rank
                if(!is.numeric.X.k & vc) {
                    ssr.rank <- ssr.rank - 1
                }
            } else if(ncol.X == 1) {
                ## With only one predictor, restricted model is
                ## unconditional mean
                ssr <- sum((y-mean(y))^2)
                ssr.rank <- 1
            } else if(ncol.X>1 & Est.ssu$num.x == 1 & is.numeric.X.k) {
                foo <- X.res
                for(i in 1:NCOL(foo)) foo[,i] <- as.numeric(foo[,i])
                ## Only one numeric predictor, rest must be factors, compute multivariate mean
                z.unique <- uniquecombs(as.matrix(foo))
                rm(foo)
                ind <-  attr(z.unique,"index")
                ind.vals <-  unique(ind)
                nrow.z.unique <- nrow(z.unique)
                mv.mean <- numeric(length=nrow.X)
                for(i in 1:nrow.z.unique) {
                    zz <- ind == ind.vals[i]
                    mv.mean[zz] <- mean(y[zz])
                }     
                ssr <- sum((y-mv.mean)^2)
                ssr.rank <- 1
            }

            F.stat[k] <- (nrow.X-ssu.rank)*(ssr-ssu)/((ssu.rank-ssr.rank)*ssu)

            for(b in 1:B) {
                if(verbose) cat(paste("\rAnova for predictor ",k," of ",ncol.X," (bootstrap replication ",b," of ",B,")",sep=""))
                ## Residual bootstrap from the null model, use
                ## original model configuration with bootstrap y

                if(ncol.X>1 & (Est.ssu$num.x>1 | (Est.ssu$num.x==1 & !is.numeric.X.k))) {
                    y.boot <- Est.ssr$fitted.values + sample(c(y-Est.ssr$fitted.values)*sqrt(nrow.X/(nrow.X-ssr.rank)),size=nrow.X,replace=TRUE)
                }  else if(ncol.X == 1) {
                    y.boot <- mean(y) + sample(y-mean(y),replace=TRUE)
                }  else if(ncol.X>1 & Est.ssu$num.x == 1 & is.numeric.X.k) {
                    y.boot <- mv.mean + sample(y-mv.mean,replace=TRUE)                    
                }

                Est.ssu.boot <- lm.ma.Est(y=y.boot,
                                          X=X,
                                          X.eval=NULL,
                                          basis=basis,
                                          compute.deriv=FALSE,
                                          deriv.order=deriv.order,
                                          degree.min=degree.min,
                                          deriv.index=deriv.index,
                                          degree.max=degree.max,
                                          lambda=lambda,
                                          segments.max=segments.max,
                                          knots=knots,
                                          S=S,
                                          method=method,
                                          ma.weights=Est.ssu$ma.weights,
                                          basis.vec=Est.ssu$basis.vec,
                                          rank.vec=Est.ssu$rank.vec,
                                          K.mat=Est.ssu$DS,
                                          weights=weights,
                                          vc=vc,
                                          verbose=FALSE,
                                          tol=tol,
                                          ...)

                ssu.boot <- sum((y.boot-Est.ssu.boot$fitted.values)^2)  
                
                if(ncol.X>1 & (Est.ssu$num.x>1 | (Est.ssu$num.x==1 & !is.numeric.X.k))) {

                    Est.ssr.boot <- lm.ma.Est(y=y.boot,
                                              X=X.res,
                                              X.eval=NULL,
                                              basis=basis,
                                              compute.deriv=FALSE,
                                              deriv.order=deriv.order,
                                              degree.min=degree.min,
                                              deriv.index=deriv.index,
                                              degree.max=degree.max,
                                              lambda=lambda,
                                              segments.max=segments.max,
                                              knots=knots,
                                              S=S,
                                              method=method,
                                              ma.weights=Est.ssu$ma.weights,
                                              basis.vec=Est.ssu$basis.vec,
                                              rank.vec=Est.ssu$rank.vec,
                                              K.mat=if(is.numeric.X.k){Est.ssu$DS[,c(-k,-(k+Est.ssu$num.x)),drop=FALSE]}else{Est.ssu$DS},
                                              weights=weights,
                                              vc=vc,
                                              verbose=FALSE,
                                              tol=tol,
                                              ...)

                    ssr.boot <- sum((y.boot-Est.ssr.boot$fitted.values)^2)         

                } else if(ncol.X == 1) {
                    ssr.boot <- sum((y.boot-mean(y.boot))^2)
                }   else if(ncol.X>1 & Est.ssu$num.x == 1 & is.numeric.X.k) {
                    for(i in 1:nrow.z.unique) {
                        zz <- ind == ind.vals[i]
                        mv.mean[zz] <- mean(y.boot[zz])
                    }     
                    ssr.boot <- sum((y.boot-mv.mean)^2)
                }

                F.boot[b] <- (nrow.X-ssu.rank)*(ssr.boot-ssu.boot)/((ssu.rank-ssr.rank)*ssu.boot)


            }

            P.vec[k] <- mean(ifelse(F.boot>F.stat[k],1,0))
            
            if(verbose) cat("\r                                                                               ")
            if(verbose) cat("\r")

        }

        Est$P.vec <- P.vec
        Est$F.stat <- F.stat
    }
    
    Est$compute.anova <- compute.anova
    Est$bootstrap.ci <- bootstrap.ci
    
    if(verbose) cat("\r                                                                               ")
    if(verbose) cat("\r")

    class(Est) <- "lm.ma"
    return(Est)

}

print.lm.ma <- function(x,
                        ...) {
    cat("Call:\n")
    print(x$call)

}

## Method for extracting vector of derivatives.

coef.lm.ma <- function(object,
                       ...) {

    if(!is.null(object$deriv)) {
        object$deriv
    } else {
        cat("\nlm.ma(...) was called with option compute.deriv=FALSE\n")
    }

}

## Method for summary.

summary.lm.ma <- function(object,
                          ...) {

    cat("Call:\n")
    print(object$call)
    cat("\nModel Averaging Linear Regression",sep="")
    cat(paste(ifelse(object$vc, " (Varying Coefficient Specification)"," (Additive Dummy Specification)"),sep=""))
    cat(paste("\nModel average criterion: ", ifelse(object$method=="jma","Jackknife (Hansen and Racine (2013))","Mallows  (Hansen (2007))"), sep=""))
    cat(paste("\nMinimum degree: ", object$degree.min, sep=""))  
    cat(paste("\nMaximum degree: ", object$degree.max, sep=""))
    cat(paste("\nBasis: ", object$basis, sep=""))  
    cat(paste("\nNumber of observations: ", object$nobs, sep=""))
    cat(paste("\nEquivalent number of parameters: ", formatC(object$ma.model.rank,format="f",digits=2), sep=""))
    cat(paste("\nResidual standard error: ", format(sqrt(sum(object$residuals^2)/(object$nobs-sum(object$rank.vec*object$ma.weights))),digits=4),
              " on ", formatC(object$nobs-sum(object$rank.vec*object$ma.weights),format="f",digits=2)," degrees of freedom",sep=""))
    cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4), sep=""))
    
    cat("\n\nNon-zero model average weights: ")
    cat(formatC(object$ma.weights[object$ma.weights>1e-05]/sum(object$ma.weights[object$ma.weights>1e-05]),format="f",digits=5))
    cat("\nNon-zero weight model ranks: ")
    cat(object$rank.vec[object$ma.weights>1e-05])
    if(object$basis=="auto") {
        cat("\nNon-zero weight model bases: ")
        cat(object$basis.vec[object$ma.weights>1e-05])    
    }
    if(object$compute.anova) {
        
        reject <- rep('', length(object$P.vec))
        reject[a <- (object$P.vec < 0.1)] <- '.'
        reject[a <- (object$P.vec < 0.05)] <- '*'
        reject[a <- (object$P.vec < 0.01)] <- '**'
        reject[a <- (object$P.vec < 0.001)] <- '***'
      
        maxNameLen <- max(nc <- nchar(nm <- names(object$X)))
        maxPvalLen <- max(ncp <- nchar(format.pval(object$P.vec)))
        maxrejLen <- max(ncr <- nchar(reject))

        cat("\n\nNonparametric significance test(s)\n")
        cat("P Value(s):", paste("\n", nm," ",
                                 blank(maxNameLen-nc),
                                 format.pval(object$P.vec),
                                 blank(maxPvalLen-ncp),
                                 " ", reject,
                                 blank(maxrejLen-ncr),
                                 " [F = ",
                                 formatC(object$F.stat,format="f",digits=3),
                                 "]",
                                 sep=''))
        
        cat("\n---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
        
    }
    cat("\n\n")
    
}

## Method for predicting given a new data frame.

predict.lm.ma <- function(object,
                          newdata=NULL,
                          ...) {
    
    if(is.null(newdata)) {
        return(fitted(object))
    } else{
        
        tt <- terms(object)
        exdat <- model.frame(delete.response(tt),newdata,xlev=object$xlevels)
        Est <- lm.ma.default(y=object$y,
                             X=object$X,
                             X.eval=exdat,
                             basis=object$basis,
                             compute.deriv=object$compute.deriv,
                             deriv.index=object$deriv.index,
                             deriv.order=object$deriv.order,
                             degree.max=object$degree.max,
                             lambda=object$lambda,
                             segments.max=object$segments.max,
                             knots=object$knots,
                             S=object$S,
                             method=object$method,
                             ma.weights=object$ma.weights,
                             basis.vec=object$basis.vec,
                             rank.vec=object$rank.vec,
                             weights=object$weights,
                             vc=object$vc,
                             verbose=object$verbose,
                             tol=object$tol,
                             ...)
        
        if(object$bootstrap.ci | object$compute.deriv) {
            return(list(fit=Est$fitted.values,
                        deriv=Est$deriv,
                        fit.low=Est$fitted.ci.l,
                        fit.up=Est$fitted.ci.u,
                        deriv.low=Est$deriv.ci.l,
                        deriv.up=Est$deriv.ci.u))
        } else {
            return(Est$fitted.values)    
        }
        
    }
    
}

## Method for plotting.

plot.lm.ma <- function(x,
                       B=99,
                       plot.deriv=FALSE,
                       plot.ci=FALSE,
                       plot.data=FALSE,
                       plot.num.eval=100,
                       plot.xtrim=0.005,
                       ...) {
    
    x$verbose <- FALSE
    x$bootstrap.ci <- plot.ci

    yname <- all.vars(x$call)[1]
    xznames <- names(x$X)

    ncol.X <- NCOL(x$X)
    
    if(!plot.deriv) {
        xeval.median <- data.frame(matrix(NA,plot.num.eval,ncol.X))
        names(xeval.median) <- xznames
        is.numeric.X <- logical(ncol.X)
        for(i in 1:ncol.X) {
            is.numeric.X[i] <- is.numeric(x$X[,i])
            xeval.median[,i] <- uocquantile(x$X[,i],prob=0.5)
        }
        if(ncol.X > 1) par(mfrow=c(2,ifelse(ncol.X %%2 == 0, ncol.X/2, (ncol.X+1)/2)))
        for(i in 1:ncol.X) {
            cat(paste("\rGenerating object ",i," of ",ncol.X," to plot...",sep=""))
            xeval <- xeval.median
            xeval[,i] <- seq(uocquantile(x$X[,i],plot.xtrim),
                             uocquantile(x$X[,i],1-plot.xtrim),
                             length=plot.num.eval)
            x$compute.deriv <- FALSE
            if(plot.data) {
                plot(xeval[,i],x$y,
                     ylab=yname,
                     xlab=xznames[i],
                     cex=0.1,
                     col="grey",
                     ...)
                foo <- predict(x,newdata=xeval,bootstrap.ci=plot.ci,B=B)    
                if(!is.list(foo)) suppressWarnings(foo$fit <- foo)
                if(is.numeric(x$X[,i])) {
                    lines(xeval[,i],foo$fit,col=1)
                } else {
                    points(xeval[,i],foo$fit,bg=1,col=1,pch=21)
                }
                if(plot.ci) {
                    if(is.numeric(x$X[,i])) {
                        lines(xeval[,i],foo$fit.low,col=2,lty=2)
                        lines(xeval[,i],foo$fit.up,col=2,lty=2)
                    } else {
                        points(xeval[,i],foo$fit.low,bg=2,col=2,pch=21)
                        points(xeval[,i],foo$fit.up,bg=2,col=2,pch=21)
                    }
                }
            } else {
                foo <- predict(x,newdata=xeval,bootstrap.ci=plot.ci,B=B)
                if(!is.list(foo)) suppressWarnings(foo$fit <- foo)
                if(!plot.ci) {
                    plot(xeval[,i],foo$fit,
                         ylab=yname,
                         xlab=xznames[i],
                         type=if(is.numeric(x$X[,i])){"l"}else{"p"},
                         ...)
                } else {
                    ylim <- range(c(foo$fit.low,foo$fit.up))
                    plot(xeval[,i],foo$fit,
                         ylab=yname,
                         xlab=xznames[i],
                         type=if(is.numeric(x$X[,i])){"l"}else{"p"},
                         ylim=ylim,
                         ...)
                    if(is.numeric(x$X[,i])) {
                        lines(xeval[,i],foo$fit.low,col=2,lty=2)
                        lines(xeval[,i],foo$fit.up,col=2,lty=2)
                    } else {
                        points(xeval[,i],foo$fit.low,bg=2,col=2,pch=21)
                        points(xeval[,i],foo$fit.up,bg=2,col=2,pch=21)                       
                    }
                }
            }
        }
        
        if(ncol.X > 1) par(mfrow=c(1,1))
        
    } else {
        
        xeval.median <- data.frame(matrix(NA,plot.num.eval,ncol.X))
        names(xeval.median) <- xznames
        is.numeric.X <- logical(ncol.X)
        for(i in 1:ncol.X) {
            is.numeric.X[i] <- is.numeric(x$X[,i])
            xeval.median[,i] <- uocquantile(x$X[,i],prob=0.5)
        }
        if(ncol.X > 1) par(mfrow=c(2,ifelse(ncol.X %%2 == 0, ncol.X/2, (ncol.X+1)/2)))
        for(i in 1:ncol.X) {
            cat(paste("\rGenerating object ",i," of ",ncol.X," to plot...",sep=""))
            xeval <- xeval.median
            xeval[,i] <- seq(uocquantile(x$X[,i],plot.xtrim),
                             uocquantile(x$X[,i],1-plot.xtrim),
                             length=plot.num.eval)
            x$compute.deriv <- TRUE
            x$deriv.index <- i
 
            foo <- predict(x,newdata=xeval,bootstrap.ci=plot.ci,B=B)
            if(!plot.ci) {
                plot(xeval[,i],foo$deriv[,1],
                     ylab=if(is.numeric(x$X[,i])){paste("d ",yname," / d ",xznames[i],sep="")}else{paste("Delta ",yname,sep="")},
                     xlab=xznames[i],
                     type=if(is.numeric(x$X[,i])){"l"}else{"p"},
                     ...)
            } else {
                ylim <- range(c(foo$deriv.low[,1],foo$deriv.up[,1]))
                plot(xeval[,i],foo$deriv[,1],
                     ylab=if(is.numeric(x$X[,i])){paste("d ",yname," / d ",xznames[i],sep="")}else{paste("Delta ",yname,sep="")},
                     xlab=xznames[i],
                     type=if(is.numeric(x$X[,i])){"l"}else{"p"},
                     ylim=ylim,
                     ...)
                if(is.numeric(x$X[,i])) {
                    lines(xeval[,i],foo$deriv.low[,1],col=2,lty=2)
                    lines(xeval[,i],foo$deriv.up[,1],col=2,lty=2)
                } else {
                    points(xeval[,i],foo$deriv.low[,1],bg=2,col=2,pch=21)
                    points(xeval[,i],foo$deriv.up[,1],bg=2,col=2,pch=21)                       
                }
            }
        }
        
        if(ncol.X > 1) par(mfrow=c(1,1))        

    }
        
    cat("\r                                                     ")
    cat("\r")
    
}

## The workhorse function.

lm.ma.Est <- function(y=NULL,
                      X=NULL,
                      X.eval=NULL,
                      basis.vec=NULL,
                      basis=c("tensor","taylor","additive","auto"),
                      compute.deriv=FALSE,
                      degree.max=NULL,
                      degree.min=1,
                      deriv.index=NULL,
                      deriv.order=1,
                      K.mat=NULL,
                      knots=FALSE,
                      lambda=1e-02,
                      ma.weights=NULL,
                      method=c("jma","mma"),
                      rank.vec=NULL,
                      S=2,
                      segments.max=3,
                      tol=1e-12,
                      vc=TRUE,
                      verbose=TRUE,
                      weights=NULL,                      
                      ...) {

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
    ncol.X <- NCOL(X)
    numeric.logical <- xztmp$numeric.logical
    xzindex <- 1:ncol.X
    xzindex[numeric.logical] <- 1:num.x
    if(!is.null(num.z)) xzindex[!numeric.logical] <- 1:num.z
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
        xeval <- x
        zeval <- z
    }

    ## Construct base levels for computing differences for factors,
    ## only one row but keep as a data frame

    if(!is.null(num.z)) {
        zeval.base <- zeval[1,,drop=FALSE]
        if(vc) {
            for(i in 1:num.z) zeval.base[1,i] <- min(z[,i])
        } else {
            for(i in 1:num.z) zeval.base[1,i] <- levels(z[,i])[1]
        }
    }

    rm(xztmp)

    nrow.X <- NROW(X)
    if(!is.null(X.eval)) nrow.Xeval <- NROW(X.eval)
    
    ## If there is only one numeric predictor, use additive basis
    ## (waste to use auto in this case as all bases coincide)
    
    if(num.x == 1 & (basis != "additive")) {
        basis <- "additive"
    }

    if(vc & !is.null(num.z)) {
        z.unique <- uniquecombs(as.matrix(z))
        num.z <- ncol(z.unique)
        ind <-  attr(z.unique,"index")
        ind.vals <-  unique(ind)
        nrow.z.unique <- nrow(z.unique)
        zeval.unique <- uniquecombs(as.matrix(zeval))
        num.zeval <- ncol(zeval.unique)
        ind.zeval <-  attr(zeval.unique,"index")
        ind.zeval.vals <-  unique(ind.zeval)
        nrow.zeval.unique <- nrow(zeval.unique)
        num.eval <- nrow(zeval)
    }

    if(is.null(z)) {
        include <- NULL
    } else {
        include <- rep(1,num.z)
        lambda.vec <- rep(lambda,num.z)
    }

    if(is.null(degree.max)) {
        degree.max <- max(2,ceiling(S*log(nrow.X)/num.x))
    }
    
    P <- degree.max
    ma.weights.orig <- ma.weights
    basis.vec.orig <- basis.vec
    rank.vec.orig <- rank.vec

    if(is.null(K.mat)) {
        if(knots) {
            K.mat <- matrix.combn(K.vec1=degree.min:degree.max,K.vec2=1:segments.max,num.x=num.x)
        } else {
            K.mat <- matrix.combn(K.vec1=degree.min:degree.max,num.x=num.x)
            K.mat <- cbind(K.mat[,1:num.x],matrix(1,nrow(K.mat),num.x,byrow=TRUE))
        }
    }
    
    K.mat.orig <- K.mat
    
    if(basis=="auto" & is.null(basis.vec) & is.null(ma.weights)) {
        basis.vec <- character()
    } else if(basis=="auto" & !is.null(basis.vec) & !is.null(ma.weights)) {
        basis.vec <- basis.vec[ma.weights>1e-05]
    }  else if(basis!="auto") {
        basis.vec <- rep(basis,nrow(K.mat))
    }
    if(!is.null(ma.weights)) {
        rank.vec <- rank.vec[ma.weights>1e-05]
        K.mat <- K.mat[ma.weights>1e-05,,drop=FALSE]
        ma.weights <- ma.weights[ma.weights>1e-05]/sum(ma.weights[ma.weights>1e-05])
    }

    P.num <- NROW(K.mat)

    if(is.null(deriv.index)) deriv.index <- 1:ncol.X
    num.deriv <- length(deriv.index)

    if(compute.deriv) {
        deriv.mat <- array(NA,c(if(is.null(X.eval)){nrow.X}else{nrow.Xeval},P.num,num.deriv))
        deriv <- matrix(NA,if(is.null(X.eval)){nrow.X}else{nrow.Xeval},num.deriv)
        colnames(deriv) <- names(X)[deriv.index]
    } else {
        deriv <- NULL
    }

    if(is.null(rank.vec)) rank.vec <- numeric(length=P.num)
    sigsq <- numeric(length=P.num)
    ma.mat <- matrix(NA,nrow.X,P.num)
    fitted.mat <- matrix(NA,if(is.null(X.eval)){nrow.X}else{nrow.Xeval},P.num)

    for(p in P.num:1) {

        DS <- cbind(K.mat[p,1:num.x],K.mat[p,(num.x+1):(2*num.x)])   

        if(verbose) cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,")...",sep=""))

        if(is.null(ma.weights)) {

            if(vc & !is.null(num.z)) {

                if(basis=="auto") {
                    cv.min <- Inf
                    for(b.basis in c("tensor","taylor","additive")) {
                        
                        fit.spline <- numeric(length=nrow.X)
                        htt <- numeric(length=nrow.X)
                        basis.singular <- logical(length=nrow.z.unique)
                        for(i in 1:nrow.z.unique) {
                            zz <- ind == ind.vals[i]
                            L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                            if(!is.null(weights)) L <- weights*L
                            P <- prod.spline(x=x,K=DS,knots="quantiles",basis=b.basis)
                            if(!is.fullrank(P)) basis.singular[i] <- TRUE
                            if(b.basis=="additive" || b.basis=="taylor") {
                                model.z.unique <- lm(y~P,weights=L)
                            } else {
                                model.z.unique <- lm(y~P-1,weights=L)
                            }
                            htt[zz] <- hatvalues(model.z.unique)[zz]
                            P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=b.basis))
                            fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                        }
                        
                        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                        cv.val <- mean((y - fit.spline)^2/(1-htt)^2)
                        
                        if(cv.val < cv.min) {
                            cv.min <- cv.val
                            fit.spline.min <- fit.spline
                            htt.min <- htt
                            rank.min <- model.z.unique$rank
                            basis.singular.min <- any(basis.singular==TRUE)
                            basis.vec[p] <- b.basis
                        }
 
                    }

                    if(basis.singular.min & verbose) warning("Dimension basis is ill-conditioned - reduce degree.max")

                    fit.spline <- fit.spline.min
                    htt.min <- htt
                    model.z.unique$rank <- rank.min

                } else {

                    fit.spline <- numeric(length=nrow.X)
                    htt <- numeric(length=nrow.X)
                    basis.singular <- logical(length=nrow.z.unique)
                    for(i in 1:nrow.z.unique) {
                        zz <- ind == ind.vals[i]
                        L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                        if(!is.null(weights)) L <- weights*L
                        P <- prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                        if(!is.fullrank(P)) basis.singular[i] <- TRUE
                        if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                            model.z.unique <- lm(y~P,weights=L)
                        } else {
                            model.z.unique <- lm(y~P-1,weights=L)
                        }
                        htt[zz] <- hatvalues(model.z.unique)[zz]
                        P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                        fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                    }
                    if(any(basis.singular==TRUE) & verbose) warning("Dimension basis is ill-conditioned - reduce degree.max")
                 }

                fitted.mat[,p] <- fit.spline
                    
                if(method=="mma") {
                    ma.mat[,p] <- y - fit.spline
                } else {
                    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                    ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                }
                
                rank.vec[p] <- model.z.unique$rank
                sigsq[p] <- sqrt(sum((y - fit.spline)^2)/(nrow.X-model.z.unique$rank))

            } else {

                if(basis=="auto") {
                    cv.min <- Inf
                    basis.singular <- logical(1)
                    for(b.basis in c("tensor","taylor","additive")) {
                        P <- prod.spline(x=x,z=z,K=DS,I=include,knots="quantiles",basis=b.basis)
                        if(!is.fullrank(P)) basis.singular <- TRUE
                        if(b.basis=="additive" || b.basis=="taylor") {
                            model.ma <- lm(y~P,weights=weights)
                        } else {
                            model.ma <- lm(y~P-1,weights=weights)
                        }
                        htt <- hatvalues(model.ma)
                        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                        cv.val <- mean((y - fitted(model.ma))^2/(1-htt)^2)
                        if(cv.val < cv.min) {
                            cv.min <- cv.val
                            fit.spline.min <- fitted(model.ma)
                            model.ma.min <- model.ma
                            basis.vec[p] <- b.basis
                            basis.singular.min <- basis.singular
                        }
                    }
                    
                    if(basis.singular.min & verbose) warning("Dimension basis is ill-conditioned - reduce degree.max")
                    
                    model.ma <- model.ma.min
                    fit.spline <- fitted(model.ma)

                } else {

                    P <- prod.spline(x=x,z=z,K=DS,I=include,knots="quantiles",basis=basis.vec[p])
                    if(!is.fullrank(P) & verbose) warning("Dimension basis is ill-conditioned - reduce degree.max")
                    if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                        model.ma <- lm(y~P,weights=weights)
                    } else {
                        model.ma <- lm(y~P-1,weights=weights)
                    }
                    fit.spline <- fitted(model.ma)
                }
                
                fitted.mat[,p] <- fit.spline

                if(method=="mma") {
                    ma.mat[,p] <- y - fit.spline
                } else {
                    htt <- hatvalues(model.ma)
                    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                    ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                }

                rank.vec[p] <- model.ma$rank

                sigsq[p] <- sqrt(sum(residuals(model.ma)^2)/(nrow.X-model.ma$rank))

            }
            
        }

        ## Evaluation matters - xeval and zeval set to x and z if no evaluation so use throughout

        if(!is.null(ma.weights))  {
            if(vc & !is.null(num.z)) {

                fit.spline <- numeric(length=num.eval)
                for(i in 1:nrow.zeval.unique) {
                    zz <- ind.zeval == ind.zeval.vals[i]
                    L <- suppressWarnings(prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z))
                    if(!is.null(weights)) L <- weights*L
                    P <- prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                    if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                        model.z.unique <- lm(y~P,weights=L)
                    } else {
                        model.z.unique <- lm(y~P-1,weights=L)
                    }
                    P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                    fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))

                }

                fitted.mat[,p] <- fit.spline
                rank.vec[p] <- model.z.unique$rank
                
            } else {

                P <- prod.spline(x=x,z=z,K=DS,I=include,knots="quantiles",basis=basis.vec[p])
                if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                    model.ma <- lm(y~P,weights=weights)
                } else {
                    model.ma <- lm(y~P-1,weights=weights)
                }
                
                P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p]))

                fitted.mat[,p] <- suppressWarnings(predict(model.ma,newdata=data.frame(as.matrix(P))))
                rank.vec[p] <- model.ma$rank
                
            }
            
            if(compute.deriv) {

                if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                    K.additive <- DS
                    K.additive[,2] <- ifelse(DS[,1]==0,0,DS[,2])
                    K.additive[,1] <- ifelse(DS[,1]>0,DS[,1]-1,DS[,1])
                }

                for(k in 1:num.deriv) {

                    if(vc & !is.null(num.z)) {
                        
                        if(numeric.logical[deriv.index[k]]) {

                            ## numeric

                            deriv.spline <- numeric(length(num.eval))
                            for(i in 1:nrow.zeval.unique) {
                                zz <- ind.zeval == ind.zeval.vals[i]
                                L <- suppressWarnings(prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z))
                                if(!is.null(weights)) L <- weights*L
                                P <- prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                                P.deriv <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p],deriv.index=xzindex[deriv.index[k]],deriv=deriv.order))
                                if(basis.vec[p]=="additive") {
                                    model <- lm(y~P,weights=L)
                                    dim.P.deriv <- sum(K.additive[xzindex[deriv.index[k]],])
                                    deriv.start <- ifelse(xzindex[deriv.index[k]]!=1,sum(K.additive[1:(xzindex[deriv.index[k]]-1),])+1,1)
                                    deriv.end <- deriv.start+sum(K.additive[xzindex[deriv.index[k]],])-1
                                    deriv.ind.vec <- deriv.start:deriv.end
                                    deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
                                } else if(basis.vec[p]=="tensor") {
                                    model <- lm(y~P-1,weights=L)
                                    deriv.spline[zz] <- P.deriv%*%coef(model)
                                } else if(basis.vec[p]=="taylor") {
                                    model <- lm(y~P,weights=L)
                                    deriv.spline[zz] <- P.deriv%*%coef(model)[-1]
                                }
                            }
                            
                            model.deriv <- deriv.spline

                        } else {

                            ## factor, need base levels

                            zeval.unique.tmp <- zeval.unique
                            zeval.unique.tmp[,xzindex[deriv.index[k]]] <- zeval.base[,xzindex[deriv.index[k]]]

                            fit.spline <- numeric(length=num.eval)
                            for(i in 1:nrow.zeval.unique) {
                                zz <- ind.zeval == ind.zeval.vals[i]
                                L <- suppressWarnings(prod.kernel(Z=z,z=zeval.unique.tmp[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z))
                                if(!is.null(weights)) L <- weights*L
                                P <- prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                                if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                                    model.z.unique <- lm(y~P,weights=L)
                                } else {
                                    model.z.unique <- lm(y~P-1,weights=L)
                                }
                                P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                                
                            }
                            model.deriv <- fitted.mat[,p] - fit.spline
                        }

                    } else {

                        if(numeric.logical[deriv.index[k]]) {

                            ## numeric

                            P <- prod.spline(x=x,z=z,K=DS,I=include,knots="quantiles",basis=basis.vec[p])
                            P.deriv <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p],deriv.index=xzindex[deriv.index[k]],deriv=deriv.order))
                            dim.P.tensor <- NCOL(P)
                            deriv.ind.vec <- logical(length=NCOL(P))
                            coef.vec.model <- numeric(length=NCOL(P))
                            if(basis.vec[p]=="additive") {
                                model <- lm(y~P,weights=weights)
                                coef.vec.model <- coef(model)[-1]
                                dim.P.deriv <- sum(K.additive[xzindex[deriv.index[k]],])
                                deriv.start <- ifelse(xzindex[deriv.index[k]]!=1,sum(K.additive[1:(xzindex[deriv.index[k]]-1),])+1,1)
                                deriv.end <- deriv.start+sum(K.additive[xzindex[deriv.index[k]],])-1
                                deriv.ind.vec[deriv.start:deriv.end] <- TRUE
                            } else if(basis.vec[p]=="tensor") {
                                model <- lm(y~P-1,weights=weights)
                                coef.vec.model <- coef(model)
                                deriv.ind.vec[1:dim.P.tensor] <- TRUE
                            } else if(basis.vec[p]=="taylor") {
                                model <- lm(y~P,weights=weights)
                                coef.vec.model <- coef(model)[-1]
                                deriv.ind.vec[1:dim.P.tensor] <- TRUE
                            }
                            
                            model.deriv <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%coef.vec.model[deriv.ind.vec]

                        } else {

                            ### need base levels

                            P <- prod.spline(x=x,z=z,K=DS,I=include,knots="quantiles",basis=basis.vec[p])
                            if(basis.vec[p]=="additive" || basis.vec[p]=="taylor") {
                                model.ma <- lm(y~P,weights=weights)
                            } else {
                                model.ma <- lm(y~P-1,weights=weights)
                            }

                            zeval.base.tmp <- zeval
                            zeval.base.tmp[,xzindex[deriv.index[k]]] <- zeval.base[1,xzindex[deriv.index[k]]]
                            
                            P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include,xeval=xeval,zeval=zeval.base.tmp,knots="quantiles",basis=basis.vec[p]))

                            fit.spline <- suppressWarnings(predict(model.ma,newdata=data.frame(as.matrix(P))))
                            ## factor

                            model.deriv <- fitted.mat[,p] - fit.spline
                        }

                    }
                    
                    deriv.mat[,p,k] <- model.deriv

                }
            }

        }
            
    }
    
    
    if(is.null(ma.weights)) {

        if(verbose) cat("\r                                                    ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing model average weights...")

        ## Solve the quadratic program for the Mallows model average
        ## weights
        M <- ncol(ma.mat)
        D <- t(ma.mat)%*%ma.mat
        tol.ridge <- tol
        singular.D <- FALSE
        while(qr(D)$rank<M) {
            D <- D + diag(tol.ridge,M,M)
            tol.ridge <- tol.ridge*10
            if(verbose) warning(paste("Shrinkage factor added to D in solve.QP to ensure full rank (",tol.ridge,")",sep=""))
            singular.D <- TRUE
        }
        A <- cbind(rep(1,M),diag(1,M,M))
        b0 <- c(1,rep(0,M))
        if(method=="mma") {
            d <- -sigsq[which.max(rank.vec)]*rank.vec
        } else {
            d <- t(y)%*%ma.mat
        }        
        b <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution

        num.attempts <- 0
        while(singular.D & num.attempts < 10) {
            num.attempts <- num.attempts + 1
            ## Re-solve the quadratic program for the non-zero Mallows
            ## model average weights (trivial overhead and can only
            ## improve upon the existing weights when D is not
            ## well-conditioned)
            ma.mat.reb <- ma.mat[,b>1e-05,drop=FALSE]
            M <- ncol(ma.mat.reb)
            D <- t(ma.mat.reb)%*%ma.mat.reb
            tol.ridge <- tol
            singular.D <- FALSE
            while(qr(D)$rank<M) {
                D <- D + diag(tol.ridge,M,M)
                tol.ridge <- tol.ridge*10
                if(verbose) warning(paste("Shrinkage factor added to D in solve.QP to ensure full rank when rebalancing (",tol.ridge,")",sep=""))
                singular.D <- TRUE
            } 
            A <- cbind(rep(1,M),diag(1,M,M))
            b0 <- c(1,rep(0,M))
            rank.vec.reb <- rank.vec[b>1e-05]
            if(method=="mma") {
                d <- -sigsq[which.max(rank.vec)]*rank.vec.reb
            } else {
                d <- t(y)%*%ma.mat.reb
            }        
            b.reb <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution
            
            if(!isTRUE(all.equal(as.numeric(b[b>1e-05]),as.numeric(b.reb)))) {
                if(verbose) {
                    warning(paste("Re-running solve.QP on non-zero weight models (",length(b[b>1e-05])," initial models, ",length(b.reb[b.reb>1e-05])," rebalanced ones)",sep=""))   
                    if(!isTRUE(all.equal(b[b>1e-05],b.reb[b.reb>1e-05]))) warning(all.equal(b[b>1e-05],b.reb[b.reb>1e-05]))
                }
                b[b>1e-05] <- b.reb
            }
        }
            
    } else {
        ## For bootstrapping use weights from initial call
        b <- ma.weights
    }

    if(verbose) cat("\r                                                    ")
    if(verbose) cat("\r")
    if(verbose) cat("\rComputing fitted values...")
    fitted.values <- fitted.mat%*%b
    
    if(compute.deriv) {
        if(verbose) cat("\r                                                    ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing derivatives...")
        for(k in 1:num.deriv) {
            if(length(b)>1) {
                deriv[,k] <- deriv.mat[,,k]%*%b
            } else {
                ## b must be the scalar 1
                deriv[,k] <- deriv.mat[,,k]
            }
        }
    }

    if(verbose) {
        cat("\r                                                    ")
        cat("\r")
    }

    return(list(fitted.values=fitted.values,
                deriv=deriv,
                ma.weights=if(is.null(ma.weights)){abs(b)}else{ma.weights.orig},
                basis.vec=if(is.null(ma.weights)){basis.vec}else{basis.vec.orig},
                y=y,
                X=X,
                basis=basis,
                compute.deriv=compute.deriv,
                deriv.order=deriv.order,
                degree.min=degree.min,
                deriv.index=deriv.index,
                degree.max=degree.max,
                lambda=lambda,
                segments.max=segments.max,
                knots=knots,
                S=S,
                method=method,
                num.x=num.x,
                num.z=num.z,
                numeric.logical=numeric.logical,
                rank.vec=rank.vec,
                ma.model.rank=sum(rank.vec*abs(b)),
                nobs=nrow.X,
                DS=K.mat.orig,
                vc=vc,
                verbose=verbose,
                tol=tol,
                xnames=xnames,
                znames=znames))

}

