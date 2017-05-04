## Function call.

lm.ma <- function(...) UseMethod("lm.ma")

## Default formula interface.

lm.ma.formula <- function(formula,
                          data=list(),
                          y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          all.combinations=TRUE,
                          alpha=0.05,
                          auto.basis=c("tensor","additive"),
                          auto.reduce=TRUE,
                          B=99,
                          basis.vec=NULL,
                          basis=c("auto","tensor","taylor","additive"),
                          bootstrap.ci=FALSE,
                          compute.anova=FALSE,
                          compute.anova.index=NULL,
                          compute.deriv=FALSE,
                          compute.mean=TRUE,
                          degree.by=2,
                          degree.max=NULL,
                          degree.min=0,
                          deriv.index=NULL,
                          deriv.order=1,
                          DKL.mat=NULL,
                          eps.lambda=1e-04,
                          knots=FALSE,
                          lambda.S=2,
                          lambda.num.max=NULL,
                          ma.weights=NULL,
                          ma.weights.cutoff=1e-04,
                          max.num.candidate.models=2500,
                          method=c("jma","mma"),
                          parallel=FALSE,
                          parallel.cores=NULL,
                          rank.vec=NULL,
                          restrict.sum.ma.weights=TRUE,
                          rng.seed=42,
                          S=1,
                          segments.by=2,
                          segments.max=3,
                          segments.min=1,
                          singular.ok=TRUE,
                          trace=FALSE,
                          vc=TRUE,
                          verbose=TRUE,
                          weights=NULL,
                          ...) {
    
    if(exists(".Random.seed", .GlobalEnv)) {
        save.seed <- get(".Random.seed", .GlobalEnv)
        exists.seed = TRUE
    } else {
        exists.seed = FALSE
    }
    
    set.seed(rng.seed)
    
    ptm <- proc.time()
    
    mf <- model.frame(formula=formula, data=data)
    mt <- attr(mf, "terms")
    tydat <- model.response(mf)
    txdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    Est <- lm.ma.default(y=tydat,
                         X=txdat,
                         X.eval=NULL,
                         all.combinations=all.combinations,
                         alpha=alpha,
                         auto.basis=auto.basis,
                         auto.reduce=auto.reduce,
                         B=B,
                         basis.vec=basis.vec,
                         basis=basis,
                         bootstrap.ci=bootstrap.ci,
                         compute.anova=compute.anova,
                         compute.anova.index=compute.anova.index,
                         compute.deriv=compute.deriv,
                         compute.mean=compute.mean,
                         degree.by=degree.by,
                         degree.max=degree.max,
                         degree.min=degree.min,
                         deriv.index=deriv.index,
                         deriv.order=deriv.order,
                         DKL.mat=DKL.mat,
                         knots=knots,
                         eps.lambda=eps.lambda,
                         lambda.S=lambda.S,
                         lambda.num.max=lambda.num.max,
                         ma.weights=ma.weights,
                         ma.weights.cutoff=ma.weights.cutoff,
                         max.num.candidate.models=max.num.candidate.models,
                         method=method,
                         parallel=parallel,
                         parallel.cores=parallel.cores,
                         rank.vec=rank.vec,
                         restrict.sum.ma.weights=restrict.sum.ma.weights,
                         rng.seed=rng.seed,
                         S=S,
                         segments.by=segments.by,
                         segments.max=segments.max,
                         segments.min=segments.min,
                         singular.ok=singular.ok,
                         trace=trace,
                         vc=vc,
                         verbose=verbose,
                         weights=weights,
                         ...)
    
    Est$r.squared <- RSQfunc(tydat,Est$fitted.values)
    Est$residuals <- tydat - Est$fitted.values
    
    Est$call <- match.call()
    Est$formula <- formula
    Est$terms <- mt
    Est$xlevels <- .getXlevels(mt, mf)
    
    Est$ptm <- (proc.time() - ptm)[3]
    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    return(Est)
    
}

## Default function.

lm.ma.default <- function(y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          all.combinations=TRUE,
                          alpha=0.05,
                          auto.basis=c("tensor","additive"),
                          auto.reduce=TRUE,
                          B=99,
                          basis.vec=NULL,
                          basis=c("auto","tensor","taylor","additive"),
                          bootstrap.ci=FALSE,
                          compute.anova=FALSE,
                          compute.anova.index=NULL,
                          compute.deriv=FALSE,
                          compute.mean=TRUE,
                          degree.by=2,
                          degree.max=NULL,
                          degree.min=0,
                          deriv.index=NULL,
                          deriv.order=1,
                          DKL.mat=NULL,
                          eps.lambda=1e-04,
                          knots=FALSE,
                          lambda.S=2,
                          lambda.num.max=NULL,
                          ma.weights=NULL,
                          ma.weights.cutoff=1e-04,
                          max.num.candidate.models=2500,
                          method=c("jma","mma"),
                          parallel=FALSE,
                          parallel.cores=NULL,
                          rank.vec=NULL,
                          restrict.sum.ma.weights=TRUE,
                          rng.seed=42,
                          S=1,
                          segments.by=2,
                          segments.max=3,
                          segments.min=1,
                          singular.ok=TRUE,
                          trace=FALSE,
                          vc=TRUE,
                          verbose=TRUE,
                          weights=NULL,
                          ...) {

    basis <- match.arg(basis)
    method <- match.arg(method)

    if(!is.logical(bootstrap.ci)) stop("bootstrap.ci must be either TRUE or FALSE")
    if(!is.logical(compute.deriv)) stop("compute.deriv must be either TRUE or FALSE")
    if(!is.logical(compute.mean)) stop("compute.mean must be either TRUE or FALSE")
    if(!is.logical(knots)) stop("knots must be either TRUE or FALSE")
    if(!is.logical(parallel)) stop("parallel must be either TRUE or FALSE")
    if(!is.logical(singular.ok)) stop("singular.ok must be either TRUE or FALSE")
    if(!is.logical(trace)) stop("trace must be either TRUE or FALSE")
    if(!is.logical(vc)) stop("vc must be either TRUE or FALSE")
    if(!is.logical(verbose)) stop("verbose must be either TRUE or FALSE")
    if(!is.null(X) & !is.null(X.eval) & NCOL(X)!=NCOL(X.eval)) stop("X and X.eval must contain the same number of predictors")
    if(!is.null(compute.anova.index)) if(any(compute.anova.index < 1) | any(compute.anova.index > NCOL(X))) stop("anova indices must correspond to columns of X")
    if(!is.null(degree.max)) if(degree.max < 0) stop("degree.max cannot be negative")
    if(!is.null(degree.max)) if(degree.max < degree.min) stop("degree.mix must not exceed degree.max")
    if(!is.null(deriv.index)) if(any(deriv.index < 1) | any(deriv.index > NCOL(X))) stop("derivative indices must correspond to columns of X")
    if(!is.null(lambda.num.max)) if(lambda.num.max < 1) stop("lambda.num.max must be a positive integer")
    if(!is.null(parallel.cores)) if(parallel.cores < 1) stop("the number of cores requested must be a positive integer")
    if(degree.min < 0) stop("minimum degree (degree.min) must be non-negative")
    if(eps.lambda < 0) stop("eps.lambda must be non-negative")
    if(is.null(y) | is.null(X)) stop("You must provide data for y and X")
    if(lambda.S < 1) stop("lambda.S must be a positive integer")
    if(max.num.candidate.models<1) stop("max.num.candidate.models must be a positive integer")
    if(segments.max < 1) stop("segments.max must be a positive integer")
    if(segments.max < segments.min) stop("segments.mix must not exceed segments.max")
    if(segments.min < 1) stop("segments.min must be a positive integer")
    
    ## First obtain weights, then in subsequent call computes fits and
    ## derivatives

    ## Can this be avoided on subsequent calls?
    ## if(is.null(ma.weights)) compute, otherwise pass arguments?

    Est <- NULL
    if(is.null(ma.weights)) {
        Est <- lm.ma.Est(y=y,
                         X=X,
                         X.eval=X.eval,
                         all.combinations=all.combinations,
                         auto.basis=auto.basis,
                         auto.reduce=auto.reduce,
                         basis=basis,
                         compute.deriv=FALSE,
                         compute.mean=TRUE,
                         deriv.order=deriv.order,
                         degree.min=degree.min,
                         deriv.index=deriv.index,
                         degree.by=degree.by,
                         degree.max=degree.max,
                         eps.lambda=eps.lambda,
                         lambda.S=lambda.S,
                         lambda.num.max=lambda.num.max,
                         segments.by=segments.by,
                         segments.min=segments.min,
                         segments.max=segments.max,
                         singular.ok=singular.ok,
                         trace=trace,
                         knots=knots,
                         S=S,
                         method=method,
                         ma.weights=ma.weights,
                         ma.weights.cutoff=ma.weights.cutoff,
                         max.num.candidate.models=max.num.candidate.models,
                         basis.vec=basis.vec,
                         parallel=parallel,
                         parallel.cores=parallel.cores,
                         rank.vec=rank.vec,
                         restrict.sum.ma.weights=restrict.sum.ma.weights,
                         DKL.mat=DKL.mat,
                         weights=weights,
                         vc=vc,
                         verbose=verbose,
                         ...)

        rank.vec <- Est$rank.vec ## needed?
        DKL.mat <- Est$DKL.mat
        basis.vec <- Est$basis.vec

    }

    ## Save rank vector and matrix of degrees and segments (not
    ## computed in next call since ma.weights is passed, but needed
    ## for summary)

    if(compute.deriv | !is.null(X.eval)) {

        Est <- lm.ma.Est(y=y,
                         X=X,
                         X.eval=X.eval,
                         all.combinations=all.combinations,
                         auto.basis=auto.basis,
                         auto.reduce=auto.reduce,
                         basis=basis,
                         compute.deriv=compute.deriv,
                         compute.mean=compute.mean,
                         deriv.order=deriv.order,
                         degree.min=degree.min,
                         deriv.index=deriv.index,
                         degree.by=degree.by,
                         degree.max=degree.max,
                         eps.lambda=eps.lambda,
                         lambda.S=lambda.S,
                         lambda.num.max=lambda.num.max,
                         segments.by=segments.by,
                         segments.min=segments.min,
                         segments.max=segments.max, 
                         singular.ok=singular.ok,
                         trace=trace,
                         knots=knots,
                         S=S,
                         method=method,
                         ma.weights=if(is.null(Est)){ma.weights}else{Est$ma.weights},
                         ma.weights.cutoff=ma.weights.cutoff,
                         max.num.candidate.models=max.num.candidate.models,
                         basis.vec=if(is.null(Est)){basis.vec}else{Est$basis.vec},
                         parallel=parallel,
                         parallel.cores=parallel.cores,
                         rank.vec=if(is.null(Est)){rank.vec}else{Est$rank.vec},
                         restrict.sum.ma.weights=restrict.sum.ma.weights,
                         DKL.mat=if(is.null(Est)){DKL.mat}else{Est$DKL.mat},
                         weights=weights,
                         vc=vc,
                         verbose=verbose,
                         ...)
        
        Est$rank.vec <- rank.vec
        Est$DKL.mat <- DKL.mat
        Est$basis.vec <- basis.vec

    }

    ## Bootstrap if requested uses pre-computed ma.weights
    
    if(bootstrap.ci) {

        n <- n.eval <- NROW(X) 
        if(!is.null(X.eval)) {
            n.eval <- NROW(X.eval)
        }
        
        ncol.X <- NCOL(X)
    
        boot.mat <- matrix(NA,n.eval,B)

        if(is.null(deriv.index)) deriv.index <- 1:ncol.X
        num.deriv <- length(deriv.index)

        if(compute.deriv) boot.deriv.array <- array(NA,c(n.eval,B,num.deriv))

        if(!parallel) {

            for(b in 1:B) {
                if(verbose) cat(paste("\rBootstrap replication ",b," of ",B,sep=""))
                ii <- sample(1:n,replace=TRUE)
                out.boot <- lm.ma.Est(y=y[ii],
                                      X=X[ii,],
                                      X.eval=if(is.null(X.eval)){X}else{X.eval},
                                      all.combinations=all.combinations,
                                      auto.basis=auto.basis,
                                      auto.reduce=auto.reduce,
                                      basis=basis,
                                      compute.deriv=compute.deriv,
                                      compute.mean=compute.mean,
                                      deriv.order=deriv.order,
                                      degree.min=degree.min,
                                      deriv.index=deriv.index,
                                      degree.by=degree.by,
                                      degree.max=degree.max,
                                      eps.lambda=eps.lambda,
                                      lambda.S=lambda.S,
                                      lambda.num.max=lambda.num.max,
                                      segments.by=segments.by,
                                      segments.min=segments.min,
                                      segments.max=segments.max,
                                      singular.ok=singular.ok,
                                      trace=trace,
                                      knots=knots,
                                      S=S,
                                      method=method,
                                      ma.weights=Est$ma.weights,
                                      ma.weights.cutoff=ma.weights.cutoff,
                                      max.num.candidate.models=max.num.candidate.models,
                                      basis.vec=Est$basis.vec,
                                      parallel=Est$parallel,
                                      parallel.cores=Est$parallel.cores,
                                      rank.vec=Est$rank.vec,
                                      restrict.sum.ma.weights=Est$restrict.sum.ma.weights,
                                      DKL.mat=Est$DKL.mat,
                                      weights=weights,
                                      vc=vc,
                                      verbose=FALSE,
                                      ...)
                boot.mat[,b] <- out.boot$fitted
                if(compute.deriv) for(k in 1:num.deriv) boot.deriv.array[,b,k] <- out.boot$deriv[,k]
            }

        } else {

            ## parallel

            cl<-makeCluster(if(is.null(parallel.cores)){detectCores(logical=FALSE)}else{parallel.cores})
            registerDoParallel(cl)
            

            output <- foreach(b=1:B,.verbose=FALSE) %dopar% {
                if(verbose) cat(paste("\rBootstrap replication ",b," of ",B,sep=""))
                ii <- sample(1:n,replace=TRUE)
                out.boot <- lm.ma.Est(y=y[ii],
                                      X=X[ii,],
                                      X.eval=if(is.null(X.eval)){X}else{X.eval},
                                      all.combinations=all.combinations,
                                      auto.basis=auto.basis,
                                      auto.reduce=auto.reduce,
                                      basis=basis,
                                      compute.deriv=compute.deriv,
                                      compute.mean=compute.mean,
                                      deriv.order=deriv.order,
                                      degree.min=degree.min,
                                      deriv.index=deriv.index,
                                      degree.by=degree.by,
                                      degree.max=degree.max,
                                      eps.lambda=eps.lambda,
                                      lambda.S=lambda.S,
                                      lambda.num.max=lambda.num.max,
                                      segments.by=segments.by,
                                      segments.min=segments.min,
                                      segments.max=segments.max,
                                      singular.ok=singular.ok,
                                      trace=trace,
                                      knots=knots,
                                      S=S,
                                      method=method,
                                      ma.weights=Est$ma.weights,
                                      ma.weights.cutoff=ma.weights.cutoff,
                                      max.num.candidate.models=max.num.candidate.models,
                                      basis.vec=Est$basis.vec,
                                      parallel=FALSE,
                                      parallel.cores=Est$parallel.cores,
                                      rank.vec=Est$rank.vec,
                                      restrict.sum.ma.weights=Est$restrict.sum.ma.weights,
                                      DKL.mat=Est$DKL.mat,
                                      weights=weights,
                                      vc=vc,
                                      verbose=FALSE,
                                      ...)
                if(compute.deriv) {
                    for(k in 1:num.deriv) boot.deriv.array[,b,k] <- out.boot$deriv[,k]
                    list(boot.mat=out.boot$fitted,
                         boot.deriv.array=boot.deriv.array[,b,])
                } else {
                    list(boot.mat=out.boot$fitted)
                }
            }
            
            for(b in 1:B) {
                boot.mat[,b] <- output[[b]]$boot.mat
                if(compute.deriv) boot.deriv.array[,b,] <- output[[b]]$boot.deriv.array
            }

            stopCluster(cl)
            
        }
        
        if(verbose) cat("\r                                                                                               ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing quantiles...")

        Est$fitted.ci.l <- apply(boot.mat,1,quantile,prob=alpha/2,type=1,na.rm=TRUE)
        Est$fitted.ci.u <- apply(boot.mat,1,quantile,prob=1-alpha/2,type=1,na.rm=TRUE)
        Est$fitted.scale <- apply(boot.mat,1,mad,na.rm=TRUE)
        
        if(compute.deriv) {
            Est$deriv.ci.l <- Est$deriv.ci.u <- Est$deriv.scale <- matrix(NA,n.eval,num.deriv)
            for(k in 1:num.deriv) {
                Est$deriv.ci.l[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=alpha/2,type=1,na.rm=TRUE)
                Est$deriv.ci.u[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=1-alpha/2,type=1,na.rm=TRUE)     
                Est$deriv.scale[,k] <- apply(boot.deriv.array[,,k],1,mad,na.rm=TRUE)  
            }
        }
        if(verbose) cat("\r                                                                                               ")
        if(verbose) cat("\r")

    }

    if(compute.anova) {

        if(Est$num.x==1 & !is.null(Est$num.z) & vc) stop("The combination vc=TRUE with factors and only one numeric predictor is not feasible")

        if(is.null(compute.anova.index)) {
            compute.anova.index <- 1:NCOL(X)
            num.tests <- NCOL(X)
        } else {
            num.tests <- length(compute.anova.index)
        }

        numeric.logical <- Est$numeric.logical
        num.x <- Est$num.x
        num.z <- Est$num.z
        num.obs <- Est$num.obs
        xzindex <- 1:num.x
        xzindex[numeric.logical] <- 1:num.x
        num.pred <- num.x
        if(!is.null(num.z)) {
            num.pred <- num.x+num.z
            xzindex <- 1:num.pred
            xzindex[numeric.logical] <- 1:num.x
            xzindex[!numeric.logical] <- 1:num.z            
        }

        P.vec <- numeric(length=num.tests)
        F.stat <- numeric(length=num.tests)
        F.boot <- numeric(length=B)
        
        if(!is.null(X.eval)) {
            Est.ssu <- lm.ma.Est(y=y,
                                 X=X,
                                 X.eval=NULL,
                                 all.combinations=all.combinations,
                                 auto.basis=auto.basis,
                                 auto.reduce=auto.reduce,
                                 basis=basis,
                                 compute.deriv=FALSE,
                                 compute.mean=TRUE,
                                 deriv.order=deriv.order,
                                 degree.min=degree.min,
                                 deriv.index=deriv.index,
                                 degree.by=degree.by,
                                 degree.max=degree.max,
                                 eps.lambda=eps.lambda,
                                 lambda.S=lambda.S,
                                 lambda.num.max=lambda.num.max,
                                 segments.by=segments.by,
                                 segments.min=segments.min,
                                 segments.max=segments.max,
                                 singular.ok=singular.ok,
                                 trace=trace,
                                 knots=knots,
                                 S=S,
                                 method=method,
                                 ma.weights=Est$ma.weights,
                                 ma.weights.cutoff=ma.weights.cutoff,
                                 max.num.candidate.models=max.num.candidate.models,
                                 basis.vec=Est$basis.vec,
                                 parallel=Est$parallel,
                                 parallel.cores=Est$parallel.cores,
                                 rank.vec=Est$rank.vec,
                                 restrict.sum.ma.weights=Est$restrict.sum.ma.weights,
                                 DKL.mat=Est$DKL.mat,
                                 weights=weights,
                                 vc=vc,
                                 verbose=FALSE,
                                 ...)
        } else {
            Est.ssu <- Est
        }

        ssu <- sum((y-Est.ssu$fitted.values)^2) 
        ssu.rank <- Est.ssu$ma.model.rank

        for(k in 1:num.tests) {

            if(verbose & num.tests>1) cat(paste("\rAnova for predictor ",compute.anova.index[k]," of ",num.tests,sep=""))
            if(verbose & num.tests==1) cat(paste("\rAnova for predictor ",compute.anova.index[k],sep=""))

            is.numeric.X.k <- is.numeric(X[,compute.anova.index[k]])
            
            X.res <- X[,-compute.anova.index[k],drop=FALSE]

            
            ## With > 1 predictor, restricted model does not incorporate the kth predictor

            if(num.pred>1 & (num.x>1 | (num.x==1 & !is.numeric.X.k))) {

                Est.ssr <- lm.ma.Est(y=y,
                                     X=X.res,
                                     X.eval=NULL,
                                     all.combinations=all.combinations,
                                     auto.basis=auto.basis,
                                     auto.reduce=auto.reduce,
                                     basis=basis,
                                     compute.deriv=FALSE,
                                     compute.mean=TRUE,
                                     deriv.order=deriv.order,
                                     degree.min=degree.min,
                                     deriv.index=NULL,
                                     degree.by=degree.by,
                                     degree.max=degree.max,
                                     eps.lambda=eps.lambda,
                                     lambda.S=lambda.S,
                                     lambda.num.max=lambda.num.max,
                                     segments.by=segments.by,
                                     segments.min=segments.min,
                                     segments.max=segments.max,
                                     singular.ok=singular.ok,
                                     trace=trace,
                                     knots=knots,
                                     S=S,
                                     method=method,
                                     ma.weights=Est.ssu$ma.weights,
                                     ma.weights.cutoff=ma.weights.cutoff,
                                     max.num.candidate.models=max.num.candidate.models,
                                     basis.vec=Est.ssu$basis.vec,
                                     parallel=Est.ssu$parallel,
                                     parallel.cores=Est.ssu$parallel.cores,
                                     rank.vec=Est.ssu$rank.vec,
                                     restrict.sum.ma.weights=Est.ssu$restrict.sum.ma.weights,
                                     DKL.mat=if(is.numeric.X.k){Est.ssu$DKL.mat[,c(-xzindex[compute.anova.index[k]],-(xzindex[compute.anova.index[k]]+num.x)),drop=FALSE]}else{Est.ssu$DKL.mat},
                                     weights=weights,
                                     vc=vc,
                                     verbose=FALSE,
                                     ...)

                ssr <- sum((y-Est.ssr$fitted.values)^2) 
                ssr.rank <- Est.ssr$ma.model.rank
                if(!is.numeric.X.k & vc) {
                    ssr.rank <- ssr.rank - 1
                }
            } else if(num.pred == 1) {
                ## With only one predictor, restricted model is
                ## unconditional mean
                ssr <- sum((y-mean(y))^2)
                ssr.rank <- 1
            } else if(num.pred>1 & num.x == 1 & is.numeric.X.k) {
                foo <- X.res
                for(i in 1:NCOL(foo)) foo[,i] <- as.numeric(foo[,i])
                ## Only one numeric predictor, rest must be factors, compute multivariate mean
                z.unique <- uniquecombs(as.matrix(foo))
                rm(foo)
                ind <-  attr(z.unique,"index")
                ind.vals <-  unique(ind)
                nrow.z.unique <- nrow(z.unique)
                mv.mean <- numeric(length=num.obs)
                for(i in 1:nrow.z.unique) {
                    zz <- ind == ind.vals[i]
                    mv.mean[zz] <- mean(y[zz])
                }     
                ssr <- sum((y-mv.mean)^2)
                ssr.rank <- 1
            }

            F.stat[k] <- (num.obs-ssu.rank)*(ssr-ssu)/((ssu.rank-ssr.rank)*ssu)

            if(!parallel) {

                for(b in 1:B) {
                    if(verbose & num.tests>1) cat(paste("\rAnova for predictor ",compute.anova.index[k]," of ",num.tests," (bootstrap replication ",b," of ",B,")",sep=""))
                    if(verbose & num.tests==1) cat(paste("\rAnova for predictor ",compute.anova.index[k]," (bootstrap replication ",b," of ",B,")",sep=""))
                    ## Residual bootstrap from the null model, use
                    ## original model configuration with bootstrap y
                    
                    if(num.pred>1 & (num.x>1 | (num.x==1 & !is.numeric.X.k))) {
                        y.boot <- Est.ssr$fitted.values + sample(c(y-Est.ssr$fitted.values)*sqrt(num.obs/(num.obs-ssr.rank)),size=num.obs,replace=TRUE)
                    }  else if(num.pred == 1) {
                        y.boot <- mean(y) + sample(y-mean(y),replace=TRUE)
                    }  else if(num.pred>1 & num.x == 1 & is.numeric.X.k) {
                        y.boot <- mv.mean + sample(y-mv.mean,replace=TRUE)                    
                    }
                    
                    Est.ssu.boot <- lm.ma.Est(y=y.boot,
                                              X=X,
                                              X.eval=NULL,
                                              all.combinations=all.combinations,
                                              auto.basis=auto.basis,
                                              auto.reduce=auto.reduce,
                                              basis=basis,
                                              compute.deriv=FALSE,
                                              compute.mean=TRUE,
                                              deriv.order=deriv.order,
                                              degree.min=degree.min,
                                              deriv.index=NULL,
                                              degree.by=degree.by,
                                              degree.max=degree.max,
                                              eps.lambda=eps.lambda,
                                              lambda.S=lambda.S,
                                              lambda.num.max=lambda.num.max,
                                              segments.by=segments.by,
                                              segments.min=segments.min,
                                              segments.max=segments.max,
                                              singular.ok=singular.ok,
                                              trace=trace,
                                              knots=knots,
                                              S=S,
                                              method=method,
                                              ma.weights=ma.weights,
                                              ma.weights.cutoff=ma.weights.cutoff,
                                              max.num.candidate.models=max.num.candidate.models,
                                              basis.vec=basis.vec,
                                              parallel=parallel,
                                              parallel.cores=parallel.cores,
                                              rank.vec=rank.vec,
                                              restrict.sum.ma.weights=restrict.sum.ma.weights,
                                              DKL.mat=DKL.mat,
                                              weights=weights,
                                              vc=vc,
                                              verbose=FALSE,
                                              ...)
                    
                    ssu.boot <- sum((y.boot-Est.ssu.boot$fitted.values)^2)  
                    
                    if(num.pred>1 & (num.x>1 | (num.x==1 & !is.numeric.X.k))) {
                        
                        Est.ssr.boot <- lm.ma.Est(y=y.boot,
                                                  X=X.res,
                                                  X.eval=NULL,
                                                  all.combinations=all.combinations,
                                                  auto.basis=auto.basis,
                                                  auto.reduce=auto.reduce,
                                                  basis=basis,
                                                  compute.deriv=FALSE,
                                                  compute.mean=TRUE,
                                                  deriv.order=deriv.order,
                                                  degree.min=degree.min,
                                                  deriv.index=NULL,
                                                  degree.by=degree.by,
                                                  degree.max=degree.max,
                                                  eps.lambda=eps.lambda,
                                                  lambda.S=lambda.S,
                                                  lambda.num.max=lambda.num.max,
                                                  segments.by=segments.by,
                                                  segments.min=segments.min,
                                                  segments.max=segments.max,
                                                  singular.ok=singular.ok,
                                                  trace=trace,
                                                  knots=knots,
                                                  S=S,
                                                  method=method,
                                                  ma.weights=Est.ssu.boot$ma.weights,
                                                  ma.weights.cutoff=ma.weights.cutoff,
                                                  max.num.candidate.models=max.num.candidate.models,
                                                  basis.vec=Est.ssu.boot$basis.vec,
                                                  parallel=Est.ssu.boot$parallel,
                                                  parallel.cores=Est.ssu.boot$parallel.cores,
                                                  rank.vec=Est.ssu.boot$rank.vec,
                                                  restrict.sum.ma.weights=Est.ssu.boot$restrict.sum.ma.weights,
                                                  DKL.mat=if(is.numeric.X.k){Est.ssu.boot$DKL.mat[,c(-xzindex[compute.anova.index[k]],-(xzindex[compute.anova.index[k]]+num.x)),drop=FALSE]}else{Est.ssu.boot$DKL.mat},
                                                  weights=weights,
                                                  vc=vc,
                                                  verbose=FALSE,
                                                  ...)
                        
                        ssr.boot <- sum((y.boot-Est.ssr.boot$fitted.values)^2)         
                        
                    } else if(num.pred == 1) {
                        ssr.boot <- sum((y.boot-mean(y.boot))^2)
                    }   else if(num.pred>1 & num.x == 1 & is.numeric.X.k) {
                        for(i in 1:nrow.z.unique) {
                            zz <- ind == ind.vals[i]
                            mv.mean[zz] <- mean(y.boot[zz])
                        }     
                        ssr.boot <- sum((y.boot-mv.mean)^2)
                    }
                    
                    F.boot[b] <- (num.obs-ssu.rank)*(ssr.boot-ssu.boot)/((ssu.rank-ssr.rank)*ssu.boot)

                }

            } else {

                ## parallel

                cl<-makeCluster(if(is.null(parallel.cores)){detectCores(logical=FALSE)}else{parallel.cores})
                registerDoParallel(cl)
                
                output <- foreach(b=1:B,.verbose=FALSE) %dopar% {
                
                    if(verbose) cat(paste("\rAnova for predictor ",compute.anova.index[k]," of ",num.tests," (bootstrap replication ",b," of ",B,")",sep=""))
                    ## Residual bootstrap from the null model, use
                    ## original model configuration with bootstrap y
                    
                    if(num.pred>1 & (num.x>1 | (num.x==1 & !is.numeric.X.k))) {
                        y.boot <- Est.ssr$fitted.values + sample(c(y-Est.ssr$fitted.values)*sqrt(num.obs/(num.obs-ssr.rank)),size=num.obs,replace=TRUE)
                    }  else if(num.pred == 1) {
                        y.boot <- mean(y) + sample(y-mean(y),replace=TRUE)
                    }  else if(num.pred>1 & num.x == 1 & is.numeric.X.k) {
                        y.boot <- mv.mean + sample(y-mv.mean,replace=TRUE)                    
                    }
                    
                    Est.ssu.boot <- lm.ma.Est(y=y.boot,
                                              X=X,
                                              X.eval=NULL,
                                              all.combinations=all.combinations,
                                              auto.basis=auto.basis,
                                              auto.reduce=auto.reduce,
                                              basis=basis,
                                              compute.deriv=FALSE,
                                              compute.mean=TRUE,
                                              deriv.order=deriv.order,
                                              degree.min=degree.min,
                                              deriv.index=NULL,
                                              degree.by=degree.by,
                                              degree.max=degree.max,
                                              eps.lambda=eps.lambda,
                                              lambda.S=lambda.S,
                                              lambda.num.max=lambda.num.max,
                                              segments.by=segments.by,
                                              segments.min=segments.min,
                                              segments.max=segments.max,
                                              singular.ok=singular.ok,
                                              trace=trace,
                                              knots=knots,
                                              S=S,
                                              method=method,
                                              ma.weights=ma.weights,
                                              ma.weights.cutoff=ma.weights.cutoff,
                                              max.num.candidate.models=max.num.candidate.models,
                                              basis.vec=basis.vec,
                                              parallel=FALSE, ## override since already in parallel
                                              parallel.cores=parallel.cores,
                                              rank.vec=rank.vec,
                                              restrict.sum.ma.weights=restrict.sum.ma.weights,
                                              DKL.mat=DKL.mat,
                                              weights=weights,
                                              vc=vc,
                                              verbose=FALSE,
                                              ...)
                    
                    ssu.boot <- sum((y.boot-Est.ssu.boot$fitted.values)^2)  
                    
                    if(num.pred>1 & (num.x>1 | (num.x==1 & !is.numeric.X.k))) {
                        
                        Est.ssr.boot <- lm.ma.Est(y=y.boot,
                                                  X=X.res,
                                                  X.eval=NULL,
                                                  all.combinations=all.combinations,
                                                  auto.basis=auto.basis,
                                                  auto.reduce=auto.reduce,
                                                  basis=basis,
                                                  compute.deriv=FALSE,
                                                  compute.mean=TRUE,
                                                  deriv.order=deriv.order,
                                                  degree.min=degree.min,
                                                  deriv.index=NULL,
                                                  degree.by=degree.by,
                                                  degree.max=degree.max,
                                                  eps.lambda=eps.lambda,
                                                  lambda.S=lambda.S,
                                                  lambda.num.max=lambda.num.max,
                                                  segments.by=segments.by,
                                                  segments.min=segments.min,
                                                  segments.max=segments.max,
                                                  singular.ok=singular.ok,
                                                  trace=trace,
                                                  knots=knots,
                                                  S=S,
                                                  method=method,
                                                  ma.weights=Est.ssu.boot$ma.weights,
                                                  ma.weights.cutoff=ma.weights.cutoff,
                                                  max.num.candidate.models=max.num.candidate.models,
                                                  basis.vec=Est.ssu.boot$basis.vec,
                                                  parallel=FALSE,
                                                  parallel.cores=Est.ssu.boot$parallel.cores,
                                                  rank.vec=Est.ssu.boot$rank.vec,
                                                  restrict.sum.ma.weights=Est.ssu.boot$restrict.sum.ma.weights,
                                                  DKL.mat=if(is.numeric.X.k){Est.ssu.boot$DKL.mat[,c(-xzindex[compute.anova.index[k]],-(xzindex[compute.anova.index[k]]+num.x)),drop=FALSE]}else{Est.ssu.boot$DKL.mat},
                                                  weights=weights,
                                                  vc=vc,
                                                  verbose=FALSE,
                                                  ...)
                        
                        ssr.boot <- sum((y.boot-Est.ssr.boot$fitted.values)^2)         
                        
                    } else if(num.pred == 1) {
                        ssr.boot <- sum((y.boot-mean(y.boot))^2)
                    }   else if(num.pred>1 & num.x == 1 & is.numeric.X.k) {
                        for(i in 1:nrow.z.unique) {
                            zz <- ind == ind.vals[i]
                            mv.mean[zz] <- mean(y.boot[zz])
                        }     
                        ssr.boot <- sum((y.boot-mv.mean)^2)
                    }
                    
                    list(F.boot=(num.obs-ssu.rank)*(ssr.boot-ssu.boot)/((ssu.rank-ssr.rank)*ssu.boot))
                    
                }

                for(b in 1:B) {
                    F.boot[b] <- output[[b]]$F.boot
                }
                
                stopCluster(cl)
            }
            
            P.vec[k] <- mean(ifelse(F.boot>F.stat[k],1,0))
            
            if(verbose) cat("\r                                                                                               ")
            if(verbose) cat("\r")

        }

        Est$P.vec <- P.vec
        Est$F.stat <- F.stat
    }

    Est$compute.anova <- compute.anova
    Est$compute.anova.index <- compute.anova.index
    Est$bootstrap.ci <- bootstrap.ci
    Est$alpha <- alpha
    
    if(verbose) cat("\r                                                                                               ")
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
    if(object$knots) {
        cat(paste("\nMinimum number of interior knots: ", object$segments.min-1, sep=""))  
        cat(paste("\nMaximum number of interior knots: ", object$segments.max-1, sep=""))  
    }
    cat(paste("\nBasis: ", object$basis, sep=""))  
    cat(paste("\nNumber of observations: ", object$num.obs, sep=""))
    cat(paste("\nNumber of numeric predictors: ", object$num.x, sep=""))
    if(!is.null(object$num.z)) cat(paste("\nNumber of categorical predictors: ", object$num.z, sep=""))
    cat(paste("\nEquivalent number of parameters: ", formatC(object$ma.model.rank,format="f",digits=2), sep=""))
    cat(paste("\nResidual standard error: ", format(sqrt(sum(object$residuals^2)/(object$num.obs-sum(object$rank.vec*object$ma.weights))),digits=4),
              " on ", formatC(object$num.obs-sum(object$rank.vec*object$ma.weights),format="f",digits=2)," degrees of freedom",sep=""))
    cat(paste("\nMultiple R-squared: ", format(object$r.squared,digits=4), sep=""))
    cat(paste("\nEstimation time: ", formatC(object$ptm,digits=1,format="f")), " seconds",sep="")

    ma.weights <- object$ma.weights[object$ma.weights>object$ma.weights.cutoff]
    rank.vec <- object$rank.vec[object$ma.weights>object$ma.weights.cutoff]
    basis.vec <- object$basis.vec[object$ma.weights>object$ma.weights.cutoff]

    ma.weights <- ma.weights[order(rank.vec)]
    basis.vec <- basis.vec[order(rank.vec)]
    rank.vec <- rank.vec[order(rank.vec)]
    
    cat("\n\nNon-zero model average weights: ")
    cat(formatC(ma.weights,format="f",digits=5))
    cat("\nNon-zero weight model ranks: ")
    cat(rank.vec)
    if(object$basis=="auto") {
        cat("\nNon-zero weight model bases: ")
        cat(basis.vec)    
    }
    if(object$compute.anova) {
        
        reject <- rep('', length(object$P.vec))
        reject[a <- (object$P.vec < 0.1)] <- '.'
        reject[a <- (object$P.vec < 0.05)] <- '*'
        reject[a <- (object$P.vec < 0.01)] <- '**'
        reject[a <- (object$P.vec < 0.001)] <- '***'
      
        maxNameLen <- max(nc <- nchar(nm <- names(object$X[,object$compute.anova.index,drop=FALSE])))
        maxPvalLen <- max(ncp <- nchar(format.pval(object$P.vec)))
        maxrejLen <- max(ncr <- nchar(reject))

        if(length(object$compute.anova.index)==1){cat("\n\nNonparametric significance test\n")}else{cat("\n\nNonparametric significance tests\n")}
        cat(if(length(object$compute.anova.index)==1){"P Value:"}else{"P Values:"}, paste("\n", nm," ",
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

    if(object$auto.reduce.invoked) cat("Note: auto.reduce invoked  - see comments in Notes section (?lm.ma)\n\n")
    
}

## Method for predicting given a new data frame.

predict.lm.ma <- function(object,
                          newdata=NULL,
                          ...) {
    
    if(is.null(newdata)) {
        return(fitted(object))
    } else{

        Est <- lm.ma.default(y=object$y,
                             X=object$X,
                             X.eval=model.frame(delete.response(terms(object)),newdata,xlev=object$xlevels),
                             alpha=object$alpha,
                             all.combinations=object$all.combinations,
                             auto.basis=object$auto.basis,
                             auto.reduce=object$auto.reduce,
                             auto.reduce.invoked=object$auto.reduce.invoked,
                             basis.vec=object$basis.vec,
                             basis=object$basis,
                             eps.lambda=object$eps.lambda,
                             compute.deriv=object$compute.deriv,
                             compute.mean=object$compute.mean,
                             degree.min=object$degree.min,
                             degree.by=object$degree.by,
                             degree.max=object$degree.max,
                             deriv.index=object$deriv.index,
                             deriv.order=object$deriv.order,
                             knots=object$knots,
                             DKL.mat=object$DKL.mat,
                             lambda.S=object$lambda.S,
                             lambda.num.max=object$lambda.num.max,
                             ma.weights=object$ma.weights,
                             ma.weights.cutoff=object$ma.weights.cutoff,
                             max.num.candidate.models=object$max.num.candidate.models,
                             method=object$method,
                             parallel=object$parallel,
                             parallel.cores=object$parallel.cores,
                             rank.vec=object$rank.vec,
                             restrict.sum.ma.weights=object$restrict.sum.ma.weights,
                             rng.seed=object$rng.seed,
                             S=object$S,
                             segments.by=object$segments.by,
                             segments.min=object$segments.min,
                             segments.max=object$segments.max,
                             singular.ok=object$singular.ok,
                             trace=object$trace,
                             vc=object$vc,
                             verbose=object$verbose,
                             weights=object$weights,
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
                       plot.ci=FALSE,
                       plot.data=FALSE,
                       plot.deriv=FALSE,
                       plot.num.eval=250,
                       plot.rug=FALSE,
                       plot.xtrim=0.005,
                       ...) {

    if(!is.logical(plot.deriv)) stop("plot.deriv must be either TRUE or FALSE") 
    if(!is.logical(plot.ci)) stop("plot.ci must be either TRUE or FALSE")   
    if(!is.logical(plot.data)) stop("plot.data must be either TRUE or FALSE")   
    if(plot.num.eval<1) stop("plot.num.eval must be positive")
    if(plot.xtrim<0 | plot.xtrim>=0.5) stop("plot.xtrim must lie in [0,0.5)")

    x$verbose <- FALSE
    x$bootstrap.ci <- plot.ci
    numeric.logical <- x$numeric.logical

    yname <- all.vars(x$call)[1]
    xznames <- names(x$X)

    ncol.X <- NCOL(x$X)
    if(!plot.deriv) {
        ## Plot partial means
        xzeval.median <- data.frame(matrix(NA,plot.num.eval,ncol.X))
        names(xzeval.median) <- xznames
        for(i in 1:ncol.X) {
            xzeval.median[,i] <- uocquantile(x$X[,i],prob=0.5)
        }
        if(ncol.X > 1) par(mfrow=n2mfrow(ncol.X))
        for(i in 1:ncol.X) {
            xzeval <- xzeval.median
            if(numeric.logical[i]) {
              xzeval[,i] <- seq(uocquantile(x$X[,i],plot.xtrim),
                                uocquantile(x$X[,i],1-plot.xtrim),
                                length=plot.num.eval)
              xlim <- c(uocquantile(x$X[,i],plot.xtrim),
                        uocquantile(x$X[,i],1-plot.xtrim))
            } else {
                u <- sort(unique(x$X[,i]))
                xzeval[1:length(u),i] <- u
                xzeval <- xzeval[1:length(u),]
            }
            x$compute.deriv <- FALSE
            if(plot.data) {
                cat(paste("\rPlotting data for object ",i," of ",ncol.X,"...",sep=""))
                plot(if(numeric.logical[i]){x$X[x$X[,i] >= xlim[1] & x$X[,i] <= xlim[2] ,i]}else{x$X[,i]},if(numeric.logical[i]){x$y[x$X[,i] >= xlim[1] & x$X[,i] <= xlim[2]]}else{x$y},
                     xlim=if(numeric.logical[i]){xlim}else{NULL},
                     ylab=yname,
                     xlab=xznames[i],
                     cex=0.1,
                     col="grey",
                     ...)
                if(numeric.logical[i] & plot.rug) suppressWarnings(rug(x$X[x$X[,i] >= xlim[1] & x$X[,i] <= xlim[2] ,i]))
                cat(paste("\rGenerating object ",i," of ",ncol.X," to plot...",sep=""))                
                foo <- predict(x,newdata=xzeval,bootstrap.ci=plot.ci,B=B,...)    
                if(!is.list(foo)) suppressWarnings(foo$fit <- foo)
                if(numeric.logical[i]) {
                    lines(xzeval[,i],foo$fit,col=1)
                } else {
                    points(xzeval[,i],foo$fit,bg=1,col=1,pch=21)
                }
                if(plot.ci) {
                    if(numeric.logical[i]) {
                        lines(xzeval[,i],foo$fit.low,col=2,lty=2)
                        lines(xzeval[,i],foo$fit.up,col=2,lty=2)
                    } else {
                        points(xzeval[,i],foo$fit.low,bg=2,col=2,pch=21)
                        points(xzeval[,i],foo$fit.up,bg=2,col=2,pch=21)
                    }
                }
            } else {
                cat(paste("\rGenerating object ",i," of ",ncol.X," to plot...",sep=""))    
                foo <- predict(x,newdata=xzeval,bootstrap.ci=plot.ci,B=B,...)
                if(!is.list(foo)) suppressWarnings(foo$fit <- foo)
                if(!plot.ci) {
                    plot(xzeval[,i],foo$fit,
                         ylab=yname,
                         xlab=xznames[i],
                         type=if(numeric.logical[i]){"l"}else{"p"},
                         ...)
                } else {
                    ylim <- range(c(foo$fit.low,foo$fit.up))
                    plot(xzeval[,i],foo$fit,
                         ylab=yname,
                         xlab=xznames[i],
                         type=if(numeric.logical[i]){"l"}else{"p"},
                         ylim=ylim,
                         ...)
                    if(numeric.logical[i]) {
                        lines(xzeval[,i],foo$fit.low,col=2,lty=2)
                        lines(xzeval[,i],foo$fit.up,col=2,lty=2)
                    } else {
                        points(xzeval[,i],foo$fit.low,bg=2,col=2,pch=21)
                        points(xzeval[,i],foo$fit.up,bg=2,col=2,pch=21)                       
                    }
                }
            }
        }
        
        if(ncol.X > 1) par(mfrow=c(1,1))
        
    } else {

        ## Plot derivatives
        xzeval.median <- data.frame(matrix(NA,plot.num.eval,ncol.X))
        names(xzeval.median) <- xznames
        for(i in 1:ncol.X) {
            xzeval.median[,i] <- uocquantile(x$X[,i],prob=0.5)
        }
        if(ncol.X > 1) par(mfrow=n2mfrow(ncol.X))
        for(i in 1:ncol.X) {
            cat(paste("\rGenerating object ",i," of ",ncol.X," to plot...",sep=""))
            xzeval <- xzeval.median
            if(numeric.logical[i]) {
                xzeval[,i] <- seq(uocquantile(x$X[,i],plot.xtrim),
                                  uocquantile(x$X[,i],1-plot.xtrim),
                                  length=plot.num.eval)
            } else {
                u <- sort(unique(x$X[,i]))
                xzeval[1:length(u),i] <- u
                xzeval <- xzeval[1:length(u),]
            }

            x$compute.deriv <- TRUE
            x$deriv.index <- i

            foo <- predict(x,newdata=xzeval,bootstrap.ci=plot.ci,B=B,...)
            if(!plot.ci) {
                plot(xzeval[,i],foo$deriv[,1],
                     ylab=if(numeric.logical[i]){paste("d ",yname," / d ",xznames[i],sep="")}else{paste("Delta ",yname,sep="")},
                     xlab=xznames[i],
                     type=if(numeric.logical[i]){"l"}else{"p"},
                     ...)
                abline(h=0,lty=2,col="grey")
            } else {
                ylim <- range(c(foo$deriv.low[,1],foo$deriv.up[,1]))
                plot(xzeval[,i],foo$deriv[,1],
                     ylab=if(numeric.logical[i]){paste("d ",yname," / d ",xznames[i],sep="")}else{paste("Delta ",yname,sep="")},
                     xlab=xznames[i],
                     type=if(numeric.logical[i]){"l"}else{"p"},
                     ylim=ylim,
                     ...)
                abline(h=0,lty=2,col="grey")                
                if(numeric.logical[i]) {
                    lines(xzeval[,i],foo$deriv.low[,1],col=2,lty=2)
                    lines(xzeval[,i],foo$deriv.up[,1],col=2,lty=2)
                } else {
                    points(xzeval[,i],foo$deriv.low[,1],bg=2,col=2,pch=21)
                    points(xzeval[,i],foo$deriv.up[,1],bg=2,col=2,pch=21)                       
                }
            }
        }
        
        if(ncol.X > 1) par(mfrow=c(1,1))        

    }
        
    cat("\r                                                                                               ")
    cat("\r")
    
}

## The workhorse function.

lm.ma.Est <- function(y=NULL,
                      X=NULL,
                      X.eval=NULL,
                      all.combinations=TRUE,
                      auto.basis=c("tensor","additive"),
                      auto.reduce=TRUE,
                      basis.vec=NULL,
                      basis=c("auto","tensor","taylor","additive"),
                      eps.lambda=1e-04,
                      compute.deriv=FALSE,
                      compute.mean=TRUE,
                      degree.by=2,
                      degree.max=NULL,
                      degree.min=0,
                      deriv.index=NULL,
                      deriv.order=1,
                      DKL.mat=NULL,
                      knots=FALSE,
                      lambda.S=2,
                      lambda.num.max=NULL,
                      ma.weights=NULL,
                      ma.weights.cutoff=1e-04,
                      max.num.candidate.models=2500,
                      method=c("jma","mma"),
                      parallel=FALSE,
                      parallel.cores=NULL,
                      rank.vec=NULL,
                      restrict.sum.ma.weights=TRUE,
                      S=1,
                      segments.by=2,
                      segments.min=1,
                      segments.max=3,
                      singular.ok=TRUE,
                      trace=TRUE,
                      vc=TRUE,
                      verbose=TRUE,
                      weights=NULL,                      
                      ...) {

    if(verbose) {
        cat("\r                                                                                               ")
        cat("\r")
        cat("\rPre-processing data...")
    }

    ma.weights.orig <- ma.weights
    basis.vec.orig <- basis.vec
    rank.vec.orig <- rank.vec
    DKL.mat.orig <- DKL.mat

    ## Take the input data frame and split into factors and numeric,
    ## get certain quantities (num.z, num.x, numeric.logical etc.)

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
    num.eval.obs <- num.obs <- NROW(X)
    if(!is.null(X.eval)) {
        if(!vc) {
            xztmp <- splitFrame(as.data.frame(X.eval))
        } else {
            xztmp <- splitFrame(as.data.frame(X.eval),factor.to.numeric=TRUE)
        }
        xeval <- xztmp$x
        zeval <- xztmp$z
        num.eval.obs <- NROW(X.eval)
        ## Test for non-overlapping support of numeric predictors
        if(trace) for(i in 1:num.x) if(min(xeval[,i]) < min(x[,i]) | max(xeval[,i]) > max(x[,i])) warning(paste(xnames[i], " evaluation data extends beyond the support of training data - prediction unreliable",sep=""),immediate.=TRUE)
        ## Test for non-overlapping support of categorical predictors
        if(!is.null(num.z) & trace) for(i in 1:num.z) if(!all(as.matrix(zeval)[,i] %in% as.matrix(z)[,i])) warning(paste(znames[i], " evaluation and training data have at least one non-overlapping level - prediction unreliable",sep=""),immediate.=TRUE)

    } else {
        ## If there is no evaluation data set to x and z (wise?)
        xeval <- x
        zeval <- z
    }
    rm(xztmp)
    
    ## Construct base levels for computing differences for factors,
    ## only one row but keep as a data frame

    categories <- NULL
    if(!is.null(num.z)) {
        zeval.base <- zeval[1,,drop=FALSE]
        categories <- numeric(length=num.z)
        for(i in 1:num.z) categories[i] <- length(unique(z[,i]))
        if(basis=="additive") categories <- categories - 1
        if(vc) {
            for(i in 1:num.z) zeval.base[1,i] <- min(z[,i])
        } else {
            for(i in 1:num.z) zeval.base[1,i] <- levels(z[,i])[1]
        }
    }

    if(vc & !is.null(num.z)) {
        if(verbose) {
            cat("\r                                                                                               ")
            cat("\r")
            cat("\rGenerating uniquecombs()...")
        }
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
    }

    if(is.null(degree.max)) {
        degree.max <- max(2,ceiling(log(num.obs)-S*log(1+num.x)))
    }

    if(is.null(lambda.num.max)) {
        lambda.num.max <- max(1,ceiling(log(num.obs)-lambda.S*log(1+num.z)))
    }

    ## If there is only one numeric predictor, use additive basis
    ## (waste to use auto in this case as all bases coincide). Also,
    ## the data-driven rules for degree.max etc. are to control
    ## complexity in multivariate settings (num.x > 1). But with only
    ## 1 numerical predictor they are overkill, so if a user is only
    ## considering the univariate case allow for more default
    ## flexibility. In particular, halve degree.by (default 2) and
    ## double degree.max to allow for more serious univariate
    ## approximation capabilities. As a competitor for loess this will
    ## likely dominate.

    if(num.x == 1 & (basis != "additive")) {
        basis <- "additive"
        if(!is.null(degree.max)) degree.max <- 2*degree.max
        degree.by <- max(1,round(degree.by/2))
        if(trace) warning("only one numeric predictor presence, degree.max doubled and degree.by halved")
    }

    degree.max.orig <- degree.max

    if(num.x == 2 & is.null(num.z)) {
        degree.by <- max(1,round(degree.by/2))
        if(trace) warning("only two numeric predictors present, no categorical, degree.by halved")
    }

    if(!is.null(num.z)) {
        if(num.x == 2 & num.z <= 2) {
            degree.by <- max(1,round(degree.by/2))
            if(trace) warning("only two numeric predictors present, two or less categorical predictors present, degree.by halved")
        }
    }
        

    lambda.seq <- NULL
    if(vc) {
        l.pow <- 1
        lambda.seq <- seq(eps.lambda,1,length=lambda.num.max)**l.pow
        n.lambda.seq <- length(lambda.seq)
    }

    if(is.null(z) | vc) {
        include <- NULL
    } else {
        include <- rep(1,num.z)
    }

    degree.seq <- c(0,seq(1,degree.max,by=degree.by))
    if(degree.min != 0) degree.seq <- seq(degree.min,degree.max,by=degree.by)
    segments.seq <- segments.min
    if(knots) segments.seq <- seq(segments.min,segments.max,by=segments.by)

    auto.reduce.flag <- FALSE
    auto.reduce.invoked <- FALSE

    if(is.null(DKL.mat)) {

        auto.reduce.num.attempts <- 0
    
        while(is.null(DKL.mat) & !auto.reduce.flag & auto.reduce.num.attempts <= 100) {

            if(verbose) {
                cat("\r                                                                                               ")
                cat("\r")
                cat("\rGenerating DKL.mat...")
            }

            if(is.null(num.z)) {
                if(knots) {
                    if(all.combinations) {
                        P.num <- nrow.DKL.mat(num.x,num.z,degree.seq,segments.seq,lambda.seq)
                        if(P.num < max.num.candidate.models) DKL.mat <- matrix.combn(K.vec1=degree.seq,K.vec2=segments.seq,num.x=num.x)
                    } else {
                        DKL.mat <- matrix(NA,length(degree.seq),2*num.x)
                        for(i in 1:length(degree.seq)) {
                            DKL.mat[i,1:num.x] <- degree.seq[i]
                            DKL.mat[i,(num.x+1):(2*num.x)] <- segments.seq[1]
                        }
                        P.num <- NROW(DKL.mat)
                    }
                } else {
                    if(all.combinations) {
                        P.num <- nrow.DKL.mat(num.x,num.z,degree.seq,segments.seq,lambda.seq)
                        if(P.num < max.num.candidate.models) DKL.mat <- matrix.combn(K.vec1=degree.seq,K.vec2=1,num.x=num.x)
                    } else {
                        DKL.mat <- matrix(NA,length(degree.seq),2*num.x)
                        for(i in 1:length(degree.seq)) {
                            DKL.mat[i,1:num.x] <- degree.seq[i]
                            DKL.mat[i,(num.x+1):(2*num.x)] <- 1
                        }
                        P.num <- NROW(DKL.mat)
                    }
                }
            } else {
                if(knots & vc) {
                    if(all.combinations) {
                        P.num <- nrow.DKL.mat(num.x,num.z,degree.seq,segments.seq,lambda.seq)
                        if(P.num < max.num.candidate.models) DKL.mat <- matrix.combn(K.vec1=degree.seq,K.vec2=segments.seq,K.vec3=lambda.seq,num.x=num.x,num.z=num.z)
                    } else {
                        DKL.mat <- matrix(NA,length(degree.seq),2*num.x+num.z)
                        for(i in 1:length(degree.seq)) {
                            DKL.mat[i,1:num.x] <- degree.seq[i]
                            DKL.mat[i,(num.x+1):(2*num.x)] <- segments.seq[1]
                            DKL.mat[i,(2*num.x+1):(2*num.x+num.z)] <- lambda.seq[1]
                        }
                        P.num <- NROW(DKL.mat)
                    }
                } else if(!knots & vc) {
                    if(all.combinations) {
                        P.num <- nrow.DKL.mat(num.x,num.z,degree.seq,segments.seq,lambda.seq)
                        if(P.num < max.num.candidate.models) DKL.mat <- matrix.combn(K.vec1=degree.seq,K.vec2=1,K.vec3=lambda.seq,num.x=num.x,num.z=num.z)
                    } else {
                        DKL.mat <- matrix(NA,length(degree.seq),2*num.x+num.z)
                        for(i in 1:length(degree.seq)) {
                            DKL.mat[i,1:num.x] <- degree.seq[i]
                            DKL.mat[i,(num.x+1):(2*num.x)] <- 1
                            DKL.mat[i,(2*num.x+1):(2*num.x+num.z)] <- lambda.seq[1]
                        }
                        P.num <- NROW(DKL.mat)
                    }
                } else {
                    if(all.combinations) {
                        P.num <- nrow.DKL.mat(num.x,num.z,degree.seq,segments.seq,lambda.seq)
                        if(P.num < max.num.candidate.models) DKL.mat <- matrix.combn(K.vec1=degree.seq,K.vec2=1,num.x=num.x)
                    } else {
                        DKL.mat <- matrix(NA,length(degree.seq),2*num.x)
                        for(i in 1:length(degree.seq)) {
                            DKL.mat[i,1:num.x] <- degree.seq[i]
                            DKL.mat[i,(num.x+1):(2*num.x)] <- 1
                        }
                        P.num <- NROW(DKL.mat)
                    }
                }
            }
            
            if(!auto.reduce & P.num >= max.num.candidate.models) {
                stop(paste("number of candidate models (",P.num,") exceeds maximum (",max.num.candidate.models,") - see comments in Notes section (?lm.ma)",""))
            }
            
            if(auto.reduce & P.num >= max.num.candidate.models) {
                if(vc & (length(lambda.seq) > 2)) lambda.seq <- lambda.seq[-(length(lambda.seq)-1)]
                if(knots & (length(segments.seq) > 1)) segments.seq <- segments.seq[-length(segments.seq)]
                if(length(degree.seq) > 2 & (length(lambda.seq) <= 2)) degree.seq <- degree.seq[-length(degree.seq)]
                if(trace) warning(paste("number of candidate models (",P.num,") exceeds maximum (",max.num.candidate.models,") - see comments in Notes section (?lm.ma)",sep=""),immediate.=TRUE)
                if(trace) warning(paste("degree.seq = ",paste(degree.seq,collapse=","),sep=""),immediate.=TRUE)
                if(vc & trace) warning(paste("lambda.seq = ",paste(lambda.seq,collapse=","),""),immediate.=TRUE)
                if(knots & trace) warning(paste("segments.seq = ",paste(segments.seq,collapse=","),sep=""),immediate.=TRUE)
                if(verbose) warning("auto.reduce invoked (set trace=TRUE to see details - see comments in Notes section (?lm.ma))")
                DKL.mat <- NULL
                auto.reduce.num.attempts <- auto.reduce.num.attempts+1
                if(auto.reduce.num.attempts==90) {
                    ## It would be virtually impossible to trip
                    ## max.num.candidate.models when all models have
                    ## the same degree, lambda, segment etc. So as we
                    ## did in the univariate case, help things out and
                    ## try to get a richer set of models.
                    all.combinations <- FALSE
                    degree.by <- 1
                    degree.max <- degree.max.orig
                    degree.seq <- c(0,seq(1,degree.max,by=degree.by))
                    if(degree.min != 0) degree.seq <- seq(degree.min,degree.max,by=degree.by)
                    segments.seq <- segments.min
                    if(knots) segments.seq <- seq(segments.min,segments.max,by=segments.by)
                    if(verbose) warning("auto.reduce invoked, setting all.combinations=FALSE as a last resort - see comments in Notes section (?lm.ma)",immediate.=TRUE)
                }
                auto.reduce.invoked <- TRUE
            } else {
                auto.reduce.flag <- TRUE
            }
            
            if(auto.reduce.num.attempts==100) stop("auto.reduce failed - see comments in Notes section (?lm.ma))")
            
        }

    }

    if(!is.null(ma.weights)) {
        rank.vec <- rank.vec[ma.weights>ma.weights.cutoff]
        DKL.mat <- DKL.mat[ma.weights>ma.weights.cutoff,,drop=FALSE]
        basis.vec <- basis.vec[ma.weights>ma.weights.cutoff]
        ma.weights <- ma.weights[ma.weights>ma.weights.cutoff]/sum(ma.weights[ma.weights>ma.weights.cutoff])
    } else {
        if(is.null(basis.vec)) basis.vec <- rep(basis,NROW(DKL.mat))
    }

    P.num <- NROW(DKL.mat)
        
    basis.singular.vec <- logical(length=P.num)

    if(is.null(deriv.index)) deriv.index <- 1:ncol.X
    num.deriv <- length(deriv.index)

    if(compute.deriv) {
        deriv.mat <- array(NA,c(if(is.null(X.eval)){num.obs}else{num.eval.obs},P.num,num.deriv))
        deriv <- matrix(NA,if(is.null(X.eval)){num.obs}else{num.eval.obs},num.deriv)
        colnames(deriv) <- names(X)[deriv.index]
    } else {
        deriv.mat <- NULL
        deriv <- NULL
    }

    if(is.null(rank.vec)) rank.vec <- numeric(length=P.num)
    sigsq <- numeric(length=P.num)
    ma.mat <- matrix(NA,num.obs,P.num)
    fitted.mat <- matrix(NA,if(is.null(X.eval)){num.obs}else{num.eval.obs},P.num)

    if(is.null(ma.weights)) {

        ## No weights passed, compute mma/jma weights

        ## For computing model average weights, scale y to have unit
        ## variance and zero mean, then rescale prior to computing
        ## fitted values. Some datasets (india from quantreg.nonpar)
        ## caused issues with Dmat being singular and I suspected
        ## simple/same issue as those arising when using raw
        ## polynomials. Confirmed produces same result as previous
        ## code but more robust.

        y <- scale(y)

        if(!parallel) {

            for(p in P.num:1) {

                DS <- cbind(DKL.mat[p,1:num.x],DKL.mat[p,(num.x+1):(2*num.x)])
                include.vec <- include
                if(!is.null(num.z) & vc) {
                    lambda.vec <- DKL.mat[p,(2*num.x+1):(2*num.x+num.z)]
                }

                if(verbose) {
                    if(length(degree.seq) < 11) {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),", lambda = ",paste(formatC(lambda.seq,format="f",digits=3),collapse=","),")",sep=""))
                        }
                    } else {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,", lambda.num = ",n.lambda.seq,")",sep=""))
                        }
                    }
                }

                ## This function is a bit odd, but when called with
                ## weights this part is not evaluated, otherwise it is.
                
                if(vc & !is.null(num.z)) {
                    
                    ## Varying coefficient regression spline formulation -
                    ## must have categorical predictors

                    if(basis=="auto") {
                        cv.min <- Inf
                        fit.spline.min <- NULL
                        for(b.basis in auto.basis) {
                           
                            fit.spline <- numeric(length=num.obs)
                            htt <- numeric(length=num.obs)
                            basis.singular <- logical(length=nrow.z.unique)
                            for(i in 1:nrow.z.unique) {
                                dim.P <- dim.bs(basis=b.basis,kernel=TRUE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                                if(dim.P/num.obs < 0.95) {
                                    zz <- ind == ind.vals[i]
                                    L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                    if(!is.null(weights)) L <- weights*L
                                    P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=b.basis))
                                    if(attr(P,"relevant")) {
                                        P.eval <- suppressWarnings(prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=b.basis))
                                        if(b.basis=="additive" | b.basis=="taylor") {
                                            model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                            fit.spline[zz] <- cbind(1,P.eval)%*%coef(model.z.unique)
                                            if(model.z.unique$rank < dim.P+1) basis.singular[i] <- TRUE
                                        } else {
                                            model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                            fit.spline[zz] <- P.eval%*%coef(model.z.unique)
                                            if(model.z.unique$rank < dim.P) basis.singular[i] <- TRUE
                                        }
                                        htt[zz] <- hat(model.z.unique$qr)[zz]
                                    } else {
                                        model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                        fit.spline[zz] <- fitted(model.z.unique)[zz]
                                        htt[zz] <- hatvalues(model.z.unique)[zz]
                                    }
                                } else {
                                    basis.singular <- !logical(length=nrow.z.unique)
                                    cv.val <- Inf
                                    htt <- NULL
                                }
                            }
                            
                            if(!any(basis.singular==TRUE)) {
                                htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                                cv.val <- mean((y - fit.spline)^2/(1-htt)^2)
                            }
                            
                            if(cv.val < cv.min & !any(basis.singular==TRUE) & (dim.P/num.obs < 0.95)) {
                                cv.min <- cv.val
                                fit.spline.min <- fit.spline
                                htt.min <- htt
                                rank.min <- model.z.unique$rank
                                basis.singular.vec[p] <- any(basis.singular==TRUE)
                                basis.vec[p] <- b.basis
                            }
                            
                        }

                        if(is.null(fit.spline.min)) stop("all bases are ill-conditioned - reduce degree.max")
                        
                        fit.spline <- fit.spline.min
                        htt <- htt.min
                        model.z.unique$rank <- rank.min
                        
                    } else {

                        ## Varying coefficient specification, non-auto basis
                        fit.spline <- numeric(length=num.obs)
                        htt <- numeric(length=num.obs)
                        basis.singular <- logical(length=nrow.z.unique)
                        for(i in 1:nrow.z.unique) {
                            dim.P <- dim.bs(basis=basis.vec[p],kernel=TRUE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                            if(dim.P/num.obs < 0.95) {
                                zz <- ind == ind.vals[i]
                                L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                if(!is.null(weights)) L <- weights*L
                                P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                                if(attr(P,"relevant")) {
                                    P.eval <- suppressWarnings(prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                    if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                        model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                        fit.spline[zz] <- cbind(1,P.eval)%*%coef(model.z.unique)
                                        if(model.z.unique$rank < dim.P+1) basis.singular[i] <- TRUE
                                    } else {
                                        model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                        fit.spline[zz] <- P.eval%*%coef(model.z.unique)
                                        if(model.z.unique$rank < dim.P) basis.singular[i] <- TRUE
                                    }
                                    htt[zz] <- hat(model.z.unique$qr)[zz]
                                } else {
                                    model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                    htt[zz] <- hatvalues(model.z.unique)[zz]
                                    fit.spline[zz] <- fitted(model.z.unique)[zz]
                                }
                            } else {
                                basis.singular <- !logical(length=nrow.z.unique)
                            }
                        }
                        basis.singular.vec[p] <- any(basis.singular==TRUE)
                    }

                    if(basis.singular.vec[p]) stop("basis is ill-conditioned - reduce degree.max")
                    fitted.mat[,p] <- fit.spline
                    
                    if(method=="mma") {
                        ma.mat[,p] <- y - fit.spline
                    } else {
                        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                        ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                    }
                    
                    rank.vec[p] <- model.z.unique$rank
                    sigsq[p] <- sqrt(sum((y - fit.spline)^2)/(num.obs-model.z.unique$rank))
                    
                } else {

                    ## Regression spline formulation (no categorical
                    ## predictors)
                    
                    if(basis=="auto") {
                        cv.min <- Inf
                        fit.spline.min <- NULL
                        for(b.basis in auto.basis) {
                            basis.singular <- logical(1)
                            dim.P <- dim.bs(basis=b.basis,kernel=FALSE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                            if(dim.P/num.obs < 0.95) {
                                P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=b.basis))
                                if(attr(P,"relevant")) {
                                    if(b.basis=="additive" | b.basis=="taylor") {
                                        if(is.null(weights)) {
                                            model.ma <- lm.fit(cbind(1,P),y,singular.ok=singular.ok)
                                        } else {
                                            model.ma <- lm.wfit(cbind(1,P),y,weights,singular.ok=singular.ok)
                                        }
                                        htt <- hat(model.ma$qr)
                                        if(model.ma$rank < dim.P+1) basis.singular <- TRUE
                                    } else {
                                        if(is.null(weights)) {
                                            model.ma <- lm.fit(P,y,singular.ok=singular.ok)
                                        } else {
                                            model.ma <- lm.wfit(P,y,weights,singular.ok=singular.ok)
                                        }
                                        htt <- hat(model.ma$qr)
                                        if(model.ma$rank < dim.P) basis.singular <- TRUE
                                    }
                                } else {
                                    model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                                    htt <- hatvalues(model.ma)
                                }
                                htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                                cv.val <- mean((y - fitted(model.ma))^2/(1-htt)^2)
                            } else {
                                basis.singular <- !logical(1)
                                cv.val <- Inf
                                htt <- NULL
                            }
                            if(cv.val < cv.min  & !basis.singular & (dim.P/num.obs < 0.95)) {
                                cv.min <- cv.val
                                fit.spline.min <- fitted(model.ma)
                                model.ma.min <- model.ma
                                basis.vec[p] <- b.basis
                                basis.singular.vec[p] <- basis.singular
                            }
                                
                        }
                        
                        if(is.null(fit.spline.min)) stop("all bases are ill-conditioned - reduce degree.max")
                    
                        fit.spline <- fit.spline.min
                        model.ma <- model.ma.min
                        
                    } else {

                        dim.P <- dim.bs(basis=basis.vec[p],kernel=FALSE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                        if(dim.P/num.obs < 0.95) {
                            P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                            if(attr(P,"relevant")) {
                                if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                    if(is.null(weights)) {
                                        model.ma <- lm.fit(cbind(1,P),y,singular.ok=singular.ok)
                                    } else {
                                        model.ma <- lm.wfit(cbind(1,P),y,weights,singular.ok=singular.ok)
                                    }
                                    htt <- hat(model.ma$qr)
                                    if(model.ma$rank < dim.P+1) basis.singular.vec[p] <- TRUE
                                } else {
                                    if(is.null(weights)) {
                                        model.ma <- lm.fit(P,y,singular.ok=singular.ok)
                                    } else {
                                        model.ma <- lm.wfit(P,y,weights,singular.ok=singular.ok)
                                    }
                                    htt <- hat(model.ma$qr)
                                    if(model.ma$rank < dim.P) basis.singular.vec[p] <- TRUE
                                }
                            } else {
                                model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                                htt <- hatvalues(model.ma)
                            }
                            fit.spline <- fitted(model.ma)
                        } else {
                            basis.singular.vec[p] <- TRUE
                        }
                    }

                    if(basis.singular.vec[p]) stop("basis is ill-conditioned - reduce degree.max")
                    fitted.mat[,p] <- fit.spline
                    
                    if(method=="mma") {
                        ma.mat[,p] <- y - fit.spline
                    } else {
                        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                        ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                    }
                    
                    rank.vec[p] <- model.ma$rank
                    
                    sigsq[p] <- sqrt(sum(residuals(model.ma)^2)/(num.obs-model.ma$rank))
                    
                }
                
            }
            
        } else {
            
            if(verbose) cat("\r                                                                                               ")
            if(verbose) cat("\r")
            if(verbose) cat("\rGenerating candidate models...")

            ## Parallel
            
            cl<-makeCluster(if(is.null(parallel.cores)){detectCores(logical=FALSE)}else{parallel.cores})
            registerDoParallel(cl)

            ## Need p to be ascending in order for dopar inorder to function

		        output <- foreach(p=1:P.num,.verbose=FALSE) %dopar% {

                DS <- cbind(DKL.mat[p,1:num.x],DKL.mat[p,(num.x+1):(2*num.x)])
                include.vec <- include
                if(!is.null(num.z) & vc) {
                    lambda.vec <- DKL.mat[p,(2*num.x+1):(2*num.x+num.z)]
                }

                if(verbose) {
                    if(length(degree.seq) < 11) {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),", lambda = ",paste(formatC(lambda.seq,format="f",digits=3),collapse=","),")",sep=""))
                        }
                    } else {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,", lambda.num = ",n.lambda.seq,")",sep=""))
                        }
                    }
                }

                ## This function is a bit odd, but when called with
                ## weights this part is not evaluated, otherwise it is.
                
                if(vc & !is.null(num.z)) {
                    
                    ## Varying coefficient regression spline formulation -
                    ## must have categorical predictors

                    if(basis=="auto") {
                        cv.min <- Inf
                        fit.spline.min <- NULL
                        for(b.basis in auto.basis) {
                           
                            fit.spline <- numeric(length=num.obs)
                            htt <- numeric(length=num.obs)
                            basis.singular <- logical(length=nrow.z.unique)
                            for(i in 1:nrow.z.unique) {
                                dim.P <- dim.bs(basis=b.basis,kernel=TRUE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                                if(dim.P/num.obs < 0.95) {
                                    zz <- ind == ind.vals[i]
                                    L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                    if(!is.null(weights)) L <- weights*L
                                    P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=b.basis))
                                    if(attr(P,"relevant")) {
                                        P.eval <- suppressWarnings(prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=b.basis))
                                        if(b.basis=="additive" | b.basis=="taylor") {
                                            model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                            fit.spline[zz] <- cbind(1,P.eval)%*%coef(model.z.unique)
                                            if(model.z.unique$rank < dim.P+1) basis.singular[i] <- TRUE
                                        } else {
                                            model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                            fit.spline[zz] <- P.eval%*%coef(model.z.unique)
                                            if(model.z.unique$rank < dim.P) basis.singular[i] <- TRUE
                                        }
                                        htt[zz] <- hat(model.z.unique$qr)[zz]
                                    } else {
                                        model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                        fit.spline[zz] <- fitted(model.z.unique)[zz]
                                        htt[zz] <- hatvalues(model.z.unique)[zz]
                                    }
                                } else {
                                    basis.singular <- !logical(length=nrow.z.unique)
                                    cv.val <- Inf
                                    htt <- NULL
                                }
                            }
                            
                            if(!any(basis.singular==TRUE)) {
                                htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                                cv.val <- mean((y - fit.spline)^2/(1-htt)^2)
                            }
                            
                            if(cv.val < cv.min & !any(basis.singular==TRUE) & (dim.P/num.obs < 0.95)) {
                                cv.min <- cv.val
                                fit.spline.min <- fit.spline
                                htt.min <- htt
                                rank.min <- model.z.unique$rank
                                basis.singular.vec[p] <- any(basis.singular==TRUE)
                                basis.vec[p] <- b.basis
                            }
                            
                        }

                        if(is.null(fit.spline.min)) stop("all bases are ill-conditioned - reduce degree.max")
                        
                        fit.spline <- fit.spline.min
                        htt <- htt.min
                        model.z.unique$rank <- rank.min
                        
                    } else {

                        ## Varying coefficient specification, non-auto basis
                        fit.spline <- numeric(length=num.obs)
                        htt <- numeric(length=num.obs)
                        basis.singular <- logical(length=nrow.z.unique)
                        for(i in 1:nrow.z.unique) {
                            dim.P <- dim.bs(basis=basis.vec[p],kernel=TRUE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                            if(dim.P/num.obs < 0.95) {
                                zz <- ind == ind.vals[i]
                                L <- prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                if(!is.null(weights)) L <- weights*L
                                P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                                if(attr(P,"relevant")) {
                                    P.eval <- suppressWarnings(prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                    if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                        model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                        fit.spline[zz] <- cbind(1,P.eval)%*%coef(model.z.unique)
                                        if(model.z.unique$rank < dim.P+1) basis.singular[i] <- TRUE
                                    } else {
                                        model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                        fit.spline[zz] <- P.eval%*%coef(model.z.unique)
                                        if(model.z.unique$rank < dim.P) basis.singular[i] <- TRUE
                                    }
                                    htt[zz] <- hat(model.z.unique$qr)[zz]
                                } else {
                                    model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                    htt[zz] <- hatvalues(model.z.unique)[zz]
                                    fit.spline[zz] <- fitted(model.z.unique)[zz]
                                }
                            } else {
                                basis.singular <- !logical(length=nrow.z.unique)
                            }
                        }
                        basis.singular.vec[p] <- any(basis.singular==TRUE)
                    }

                    if(basis.singular.vec[p]) stop("basis is ill-conditioned - reduce degree.max")
                    fitted.mat[,p] <- fit.spline
                    
                    if(method=="mma") {
                        ma.mat[,p] <- y - fit.spline
                    } else {
                        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                        ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                    }
                    
                    rank.vec[p] <- model.z.unique$rank
                    sigsq[p] <- sqrt(sum((y - fit.spline)^2)/(num.obs-model.z.unique$rank))
                    
                } else {

                    ## Regression spline formulation (no categorical
                    ## predictors)
                    
                    if(basis=="auto") {
                        cv.min <- Inf
                        fit.spline.min <- NULL
                        for(b.basis in auto.basis) {
                            basis.singular <- logical(1)
                            dim.P <- dim.bs(basis=b.basis,kernel=FALSE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                            if(dim.P/num.obs < 0.95) {
                                P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=b.basis))
                                if(attr(P,"relevant")) {
                                    if(b.basis=="additive" | b.basis=="taylor") {
                                        if(is.null(weights)) {
                                            model.ma <- lm.fit(cbind(1,P),y,singular.ok=singular.ok)
                                        } else {
                                            model.ma <- lm.wfit(cbind(1,P),y,weights,singular.ok=singular.ok)
                                        }
                                        htt <- hat(model.ma$qr)
                                        if(model.ma$rank < dim.P+1) basis.singular <- TRUE
                                    } else {
                                        if(is.null(weights)) {
                                            model.ma <- lm.fit(P,y,singular.ok=singular.ok)
                                        } else {
                                            model.ma <- lm.wfit(P,y,weights,singular.ok=singular.ok)
                                        }
                                        htt <- hat(model.ma$qr)
                                        if(model.ma$rank < dim.P) basis.singular <- TRUE
                                    }
                                } else {
                                    model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                                    htt <- hatvalues(model.ma)
                                }
                                htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                                cv.val <- mean((y - fitted(model.ma))^2/(1-htt)^2)
                            } else {
                                basis.singular <- !logical(1)
                                cv.val <- Inf
                                htt <- NULL
                            }
                            if(cv.val < cv.min  & !basis.singular & (dim.P/num.obs < 0.95)) {
                                cv.min <- cv.val
                                fit.spline.min <- fitted(model.ma)
                                model.ma.min <- model.ma
                                basis.vec[p] <- b.basis
                                basis.singular.vec[p] <- basis.singular
                            }
                                
                        }
                        
                        if(is.null(fit.spline.min)) stop("all bases are ill-conditioned - reduce degree.max")
                    
                        fit.spline <- fit.spline.min
                        model.ma <- model.ma.min
                        
                    } else {

                        dim.P <- dim.bs(basis=basis.vec[p],kernel=FALSE,degree=DS[,1],segments=DS[,2],include=include,categories=categories)
                        if(dim.P/num.obs < 0.95) {
                            P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                            if(attr(P,"relevant")) {
                                if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                    if(is.null(weights)) {
                                        model.ma <- lm.fit(cbind(1,P),y,singular.ok=singular.ok)
                                    } else {
                                        model.ma <- lm.wfit(cbind(1,P),y,weights,singular.ok=singular.ok)
                                    }
                                    htt <- hat(model.ma$qr)
                                    if(model.ma$rank < dim.P+1) basis.singular.vec[p] <- TRUE
                                } else {
                                    if(is.null(weights)) {
                                        model.ma <- lm.fit(P,y,singular.ok=singular.ok)
                                    } else {
                                        model.ma <- lm.wfit(P,y,weights,singular.ok=singular.ok)
                                    }
                                    htt <- hat(model.ma$qr)
                                    if(model.ma$rank < dim.P) basis.singular.vec[p] <- TRUE
                                }
                            } else {
                                model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                                htt <- hatvalues(model.ma)
                            }
                            fit.spline <- fitted(model.ma)
                        } else {
                            basis.singular.vec[p] <- TRUE
                        }
                    }

                    if(basis.singular.vec[p]) stop("basis is ill-conditioned - reduce degree.max")
                    fitted.mat[,p] <- fit.spline
                    
                    if(method=="mma") {
                        ma.mat[,p] <- y - fit.spline
                    } else {
                        htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                        ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                    }
                    
                    rank.vec[p] <- model.ma$rank
                    
                    sigsq[p] <- sqrt(sum(residuals(model.ma)^2)/(num.obs-model.ma$rank))
                    
                }
                
		            list(fitted.mat=fitted.mat[,p],
		                 ma.mat=ma.mat[,p],
		                 rank.vec=rank.vec[p],
		                 basis.vec=basis.vec[p],
                     basis.singular.vec=basis.singular.vec[p],
		                 sigsq=sigsq[p])
		
		        }
		
            for(p in 1:P.num) {
                fitted.mat[,p] <- output[[p]]$fitted.mat
                ma.mat[,p] <- output[[p]]$ma.mat
                rank.vec[p] <- output[[p]]$rank.vec
                basis.vec[p] <- output[[p]]$basis.vec
                basis.singular.vec[p] <- output[[p]]$basis.singular.vec
                sigsq[p] <- output[[p]]$sigsq
            }

            stopCluster(cl)

        }

        if(verbose) cat("\r                                                                                               ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing model average weights...")

        ## Solve the quadratic program for the Mallows model average
        ## weights
        M <- ncol(ma.mat)
        D <- t(ma.mat)%*%ma.mat
        tol.ridge <- sqrt(.Machine$double.eps)
        singular.D <- FALSE
        while(qr(D)$rank<M) {
            D <- D + diag(tol.ridge,M,M)
            tol.ridge <- tol.ridge*10
            if(trace) warning(paste("Shrinkage factor added to D in solve.QP to ensure full rank (",tol.ridge,")",sep=""),immediate.=TRUE)
            singular.D <- TRUE
        }
        if(method=="mma") {
            d <- -sigsq[which.max(rank.vec)]*rank.vec
        } else {
            d <- t(y)%*%ma.mat
        }        
        if(restrict.sum.ma.weights) {
            A <- cbind(rep(1,M),diag(1,M,M))
            b0 <- c(1,rep(0,M))
            b <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution
        } else {
            A <- diag(1,M,M)
            b0 <- rep(0,M)
            b <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0)$solution
            ## Not constrained to sum to one but normalize ex-post as
            ## we construct weighted averages
            b <- b/sum(b)
        }

        num.attempts <- 0
        while(singular.D & num.attempts < 10) {
            num.attempts <- num.attempts + 1
            ## Re-solve the quadratic program for the non-zero Mallows
            ## model average weights (trivial overhead and can only
            ## improve upon the existing weights when D is not
            ## well-conditioned)
            ma.mat.reb <- ma.mat[,b>ma.weights.cutoff,drop=FALSE]
            M <- ncol(ma.mat.reb)
            D <- t(ma.mat.reb)%*%ma.mat.reb
            tol.ridge <- sqrt(.Machine$double.eps)
            singular.D <- FALSE
            while(qr(D)$rank<M) {
                D <- D + diag(tol.ridge,M,M)
                tol.ridge <- tol.ridge*10
                if(trace) warning(paste("Shrinkage factor added to D in solve.QP to ensure full rank when rebalancing (",tol.ridge,")",sep=""),immediate.=TRUE)
                singular.D <- TRUE
            } 
            if(method=="mma") {
                rank.vec.reb <- rank.vec[b>ma.weights.cutoff]
                d <- -sigsq[which.max(rank.vec)]*rank.vec.reb
            } else {
                d <- t(y)%*%ma.mat.reb
            }        
            if(restrict.sum.ma.weights) {
                A <- cbind(rep(1,M),diag(1,M,M))
                b0 <- c(1,rep(0,M))
                b.reb <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution
            } else {
                A <- diag(1,M,M)
                b0 <- rep(0,M)
                b.reb <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0)$solution
                ## Not constrained to sum to one but normalize ex-post
                ## as we construct weighted averages
                b.reb <- b.reb/sum(b.reb)
            }
            
            if(!isTRUE(all.equal(as.numeric(b[b>ma.weights.cutoff]),as.numeric(b.reb)))) {
                if(trace) {
                    warning(paste("Re-running solve.QP on non-zero weight models (",length(b[b>ma.weights.cutoff])," initial models, ",length(b.reb[b.reb>ma.weights.cutoff])," rebalanced ones)",sep=""),immediate.=TRUE)   
                    if(!isTRUE(all.equal(b[b>ma.weights.cutoff],b.reb[b.reb>ma.weights.cutoff]))) warning(all.equal(b[b>ma.weights.cutoff],b.reb[b.reb>ma.weights.cutoff]),immediate.=TRUE)
                }
                b[b>ma.weights.cutoff] <- b.reb
            }
        }

        ## Scale back to original units

        fitted.mat <- as.matrix(fitted.mat * attr(y, 'scaled:scale') + attr(y, 'scaled:center'))
        y <- as.numeric(y * attr(y, 'scaled:scale') + attr(y, 'scaled:center'))
            
    } else if(!is.null(ma.weights))  {

        if(!parallel) {

            ## Weights passed, copy to b for computation of fit and/or
            ## derivatives
    
            b <- ma.weights
    
            ## NOTE - if only the derivatives are needed, can skip
            ## computing the fitted values UNLESS derivatives are for
            ## categorical predictors.
    
            for(p in P.num:1) {
    
    
                DS <- cbind(DKL.mat[p,1:num.x],DKL.mat[p,(num.x+1):(2*num.x)])
                include.vec <- include
                if(!is.null(num.z) & vc) {
                    lambda.vec <- DKL.mat[p,(2*num.x+1):(2*num.x+num.z)]
                }
		
                if(verbose) {
                    if(length(degree.seq) < 11) {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),", lambda = ",paste(formatC(lambda.seq,format="f",digits=3),collapse=","),")",sep=""))
                        }
                    } else {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,", lambda.num = ",n.lambda.seq,")",sep=""))
                        }
                    }
                }

                if(compute.mean) {
    
                    ## Compute fitted values
                     
                    if(vc & !is.null(num.z)) {
                        
                        ## Varying coefficient regression spline formulation -
                        ## must have categorical predictors
                        
                        fit.spline <- numeric(length=num.eval.obs)
                        for(i in 1:nrow.zeval.unique) {
                            zz <- ind.zeval == ind.zeval.vals[i]
                            L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                            if(!is.null(weights)) L <- weights*L
                            P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                            if(attr(P,"relevant")) {
                                P.eval <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                    model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                    fit.spline[zz] <- cbind(1,P.eval)%*%coef(model.z.unique)
                                } else {
                                    model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                    fit.spline[zz] <- P.eval%*%coef(model.z.unique)
                                }
                            } else {
                                model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(1:num.eval.obs))[zz])
                            }
                            
                        }
                        
                        fitted.mat[,p] <- fit.spline
                        rank.vec[p] <- model.z.unique$rank
                        
                    } else {
                        
                        ## Regression spline formulation (no categorical
                        ## predictors)

                        P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                        if(attr(P,"relevant")) {
                            P.eval <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p]))
                            if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                if(is.null(weights)) {
                                    model.ma <- lm.fit(cbind(1,P),y,singular.ok=singular.ok)
                                } else {
                                    model.ma <- lm.wfit(cbind(1,P),y,weights,singular.ok=singular.ok)
                                }
                                fitted.mat[,p] <- cbind(1,P.eval)%*%coef(model.ma)
                            } else {
                                if(is.null(weights)) {
                                    model.ma <- lm.fit(P,y,singular.ok=singular.ok)
                                } else {
                                    model.ma <- lm.wfit(P,y,weights,singular.ok=singular.ok)
                                }
                                fitted.mat[,p] <- P.eval%*%coef(model.ma)
                            }
                        } else {
                            model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                            fitted.mat[,p] <- suppressWarnings(predict(model.ma,newdata=data.frame(1:num.eval.obs)))
                        }
                        rank.vec[p] <- model.ma$rank
                        
                    }
                    
                }
                
                ## Compute derivatives
    
                 if(compute.deriv) {
    
                    if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                        K.additive <- DS
                        K.additive[,2] <- ifelse(DS[,1]==0,0,DS[,2])
                        K.additive[,1] <- ifelse(DS[,1]>0,DS[,1]-1,DS[,1])
                    }
    
                    for(k in 1:num.deriv) {
    
                        kk <- xzindex[deriv.index[k]]
    
                        if(vc & !is.null(num.z)) {
                            
                            if(numeric.logical[deriv.index[k]]) {
    
                                ## Compute numeric derivatives for the
                                ## varying coefficient formulation
    
                                if(DS[kk,1] != 0) {
                                    deriv.spline <- numeric(length(num.eval.obs))
                                    for(i in 1:nrow.zeval.unique) {
                                        zz <- ind.zeval == ind.zeval.vals[i]
                                        L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                        if(!is.null(weights)) L <- weights*L
                                        P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                                        P.deriv <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p],deriv.index=kk,deriv=deriv.order))
                                        if(basis.vec[p]=="additive") {
                                            model <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                            dim.P.deriv <- sum(K.additive[kk,])
                                            deriv.start <- ifelse(kk!=1,sum(K.additive[1:(kk-1),])+1,1)
                                            deriv.end <- deriv.start+sum(K.additive[kk,])-1
                                            deriv.ind.vec <- deriv.start:deriv.end
                                            deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
                                        } else if(basis.vec[p]=="tensor") {
                                            model <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                            deriv.spline[zz] <- P.deriv%*%coef(model)
                                        } else if(basis.vec[p]=="taylor") {
                                            model <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                            deriv.spline[zz] <- P.deriv%*%coef(model)[-1]
                                        }
                                    }
                                } else {
                                    deriv.spline <- rep(0,num.eval.obs)
                                }
                                
                                model.deriv <- deriv.spline
    
                            } else {
    
                                ## Compute factor `derivatives' for the
                                ## regression spline formulation
                                ## (differences), require base levels
    
                                zeval.unique.tmp <- zeval.unique
                                zeval.unique.tmp[,kk] <- zeval.base[,kk]
    
                                fit.spline <- numeric(length=num.eval.obs)
                                for(i in 1:nrow.zeval.unique) {
                                    zz <- ind.zeval == ind.zeval.vals[i]
                                    L <- prod.kernel(Z=z,z=zeval.unique.tmp[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                    if(!is.null(weights)) L <- weights*L
                                    P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                                    if(attr(P,"relevant")) {
                                        if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                            model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                        } else {
                                            model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                        }
                                        P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                        fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                                    } else {
                                        model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                        fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(1:num.eval.obs))[zz])
                                    }
                                    
                                }
                                model.deriv <- fitted.mat[,p] - fit.spline
                            }
    
                        } else {
    
                            if(numeric.logical[deriv.index[k]]) {
    
                                ## Compute numeric derivatives
    
                                if(DS[kk,1] != 0) {
                                    P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                                    P.deriv <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p],deriv.index=kk,deriv=deriv.order))
                                    dim.P.tensor <- NCOL(P)
                                    deriv.ind.vec <- logical(length=NCOL(P))
                                    if(basis.vec[p]=="additive") {
                                        model <- lm(y~P,weights=weights,singular.ok=singular.ok)
                                        coef.model <- coef(model)[-1]
                                        dim.P.deriv <- sum(K.additive[kk,])
                                        deriv.start <- ifelse(kk!=1,sum(K.additive[1:(kk-1),])+1,1)
                                        deriv.end <- deriv.start+sum(K.additive[kk,])-1
                                        deriv.ind.vec[deriv.start:deriv.end] <- TRUE
                                    } else if(basis.vec[p]=="tensor") {
                                        model <- lm(y~P-1,weights=weights,singular.ok=singular.ok)
                                        coef.model <- coef(model)
                                        deriv.ind.vec[1:dim.P.tensor] <- TRUE
                                    } else if(basis.vec[p]=="taylor") {
                                        model <- lm(y~P,weights=weights,singular.ok=singular.ok)
                                        coef.model <- coef(model)[-1]
                                        deriv.ind.vec[1:dim.P.tensor] <- TRUE
                                    }
                                    model.deriv <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%coef.model[deriv.ind.vec]
                                } else {
                                    model.deriv <- rep(0, num.eval.obs)
                                }
                                
                            } else {
    
                                ## Compute factor `derivatives'
                                ## (differences), require base levels
    
                                P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                                if(attr(P,"relevant")) {
                                    if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                        model.ma <- lm(y~P,weights=weights,singular.ok=singular.ok)
                                    } else {
                                        model.ma <- lm(y~P-1,weights=weights,singular.ok=singular.ok)
                                    }
                                    
                                    zeval.base.tmp <- zeval
                                    zeval.base.tmp[,kk] <- zeval.base[1,kk]
                                    
                                    P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,xeval=xeval,zeval=zeval.base.tmp,knots="quantiles",basis=basis.vec[p]))
                                    
                                    fit.spline <- suppressWarnings(predict(model.ma,newdata=data.frame(as.matrix(P))))
                                } else {
                                    model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                                    fit.spline <- suppressWarnings(predict(model.ma,newdata=data.frame(1:num.eval.obs)))
                                }
                                model.deriv <- fitted.mat[,p] - fit.spline
                            }
    
                        }
                        
                        deriv.mat[,p,k] <- model.deriv
                        
                    }

                }
                
            }
    
        } else {

            ## parallel

            ## Weights passed, copy to b for computation of fit and/or
            ## derivatives
    
            b <- ma.weights
    
            ## NOTE - if only the derivatives are needed, can skip
            ## computing the fitted values UNLESS derivatives are for
            ## categorical predictors.

            cl<-makeCluster(if(is.null(parallel.cores)){detectCores(logical=FALSE)}else{parallel.cores})
            registerDoParallel(cl)

            ## Need p to be ascending in order for dopar inorder to function

		        output <- foreach(p=1:P.num,.verbose=FALSE) %dopar% {
		
                DS <- cbind(DKL.mat[p,1:num.x],DKL.mat[p,(num.x+1):(2*num.x)])
                include.vec <- include
                if(!is.null(num.z) & vc) {
                    lambda.vec <- DKL.mat[p,(2*num.x+1):(2*num.x+num.z)]
                }
		
                if(verbose) {
                    if(length(degree.seq) < 11) {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree = ",paste(degree.seq,collapse=","),", lambda = ",paste(formatC(lambda.seq,format="f",digits=3),collapse=","),")",sep=""))
                        }
                    } else {
                        if(is.null(num.z) | !vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,")",sep=""))
                        } else if(vc) {
                            cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,", lambda.num = ",n.lambda.seq,")",sep=""))
                        }
                    }
                }

                if(compute.mean) {
    
                    ## Compute fitted values
                     
                    if(vc & !is.null(num.z)) {
                        
                        ## Varying coefficient regression spline formulation -
                        ## must have categorical predictors
                        
                        fit.spline <- numeric(length=num.eval.obs)
                        for(i in 1:nrow.zeval.unique) {
                            zz <- ind.zeval == ind.zeval.vals[i]
                            L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                            if(!is.null(weights)) L <- weights*L
                            P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                            if(attr(P,"relevant")) {
                                if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                    model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                } else {
                                    model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                }
                                P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                            } else {
                                model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(1:num.eval.obs))[zz])
                            }
                            
                        }
                        
                        fitted.mat[,p] <- fit.spline
                        rank.vec[p] <- model.z.unique$rank
                        
                    } else {
                        
                        ## Regression spline formulation (no categorical
                        ## predictors)
                        
                        P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                        if(attr(P,"relevant")) {
                            if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                model.ma <- lm(y~P,weights=weights,singular.ok=singular.ok)
                            } else {
                                model.ma <- lm(y~P-1,weights=weights,singular.ok=singular.ok)
                            }
                            P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p]))
                            fitted.mat[,p] <- suppressWarnings(predict(model.ma,newdata=data.frame(as.matrix(P))))
                        } else {
                            model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                            fitted.mat[,p] <- suppressWarnings(predict(model.ma,newdata=data.frame(1:num.eval.obs)))
                        }
                        rank.vec[p] <- model.ma$rank
                        
                    }
                    
                }
                
                ## Compute derivatives
    
                 if(compute.deriv) {
    
                    if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                        K.additive <- DS
                        K.additive[,2] <- ifelse(DS[,1]==0,0,DS[,2])
                        K.additive[,1] <- ifelse(DS[,1]>0,DS[,1]-1,DS[,1])
                    }
    
                    for(k in 1:num.deriv) {
    
                        kk <- xzindex[deriv.index[k]]
    
                        if(vc & !is.null(num.z)) {
                            
                            if(numeric.logical[deriv.index[k]]) {
    
                                ## Compute numeric derivatives for the
                                ## varying coefficient formulation
    
                                if(DS[kk,1] != 0) {
                                    deriv.spline <- numeric(length(num.eval.obs))
                                    for(i in 1:nrow.zeval.unique) {
                                        zz <- ind.zeval == ind.zeval.vals[i]
                                        L <- prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                        if(!is.null(weights)) L <- weights*L
                                        P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                                        P.deriv <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p],deriv.index=kk,deriv=deriv.order))
                                        if(basis.vec[p]=="additive") {
                                            model <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                            dim.P.deriv <- sum(K.additive[kk,])
                                            deriv.start <- ifelse(kk!=1,sum(K.additive[1:(kk-1),])+1,1)
                                            deriv.end <- deriv.start+sum(K.additive[kk,])-1
                                            deriv.ind.vec <- deriv.start:deriv.end
                                            deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
                                        } else if(basis.vec[p]=="tensor") {
                                            model <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                            deriv.spline[zz] <- P.deriv%*%coef(model)
                                        } else if(basis.vec[p]=="taylor") {
                                            model <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                            deriv.spline[zz] <- P.deriv%*%coef(model)[-1]
                                        }
                                    }
                                } else {
                                    deriv.spline <- rep(0,num.eval.obs)
                                }
                                
                                model.deriv <- deriv.spline
    
                            } else {
    
                                ## Compute factor `derivatives' for the
                                ## regression spline formulation
                                ## (differences), require base levels
    
                                zeval.unique.tmp <- zeval.unique
                                zeval.unique.tmp[,kk] <- zeval.base[,kk]
    
                                fit.spline <- numeric(length=num.eval.obs)
                                for(i in 1:nrow.zeval.unique) {
                                    zz <- ind.zeval == ind.zeval.vals[i]
                                    L <- prod.kernel(Z=z,z=zeval.unique.tmp[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                                    if(!is.null(weights)) L <- weights*L
                                    P <- suppressWarnings(prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p]))
                                    if(attr(P,"relevant")) {
                                        if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                            model.z.unique <- lm.wfit(cbind(1,P),y,L,singular.ok=singular.ok)
                                        } else {
                                            model.z.unique <- lm.wfit(P,y,L,singular.ok=singular.ok)
                                        }
                                        P <- suppressWarnings(prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                                        fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                                    } else {
                                        model.z.unique <- lm(y~1,weights=L,singular.ok=singular.ok)
                                        fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(1:num.eval.obs))[zz])
                                    }
                                    
                                }
                                model.deriv <- fitted.mat[,p] - fit.spline
                            }
    
                        } else {
    
                            if(numeric.logical[deriv.index[k]]) {
    
                                ## Compute numeric derivatives
    
                                if(DS[kk,1] != 0) {
                                    P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                                    P.deriv <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p],deriv.index=kk,deriv=deriv.order))
                                    dim.P.tensor <- NCOL(P)
                                    deriv.ind.vec <- logical(length=NCOL(P))
                                    if(basis.vec[p]=="additive") {
                                        model <- lm(y~P,weights=weights,singular.ok=singular.ok)
                                        coef.model <- coef(model)[-1]
                                        dim.P.deriv <- sum(K.additive[kk,])
                                        deriv.start <- ifelse(kk!=1,sum(K.additive[1:(kk-1),])+1,1)
                                        deriv.end <- deriv.start+sum(K.additive[kk,])-1
                                        deriv.ind.vec[deriv.start:deriv.end] <- TRUE
                                    } else if(basis.vec[p]=="tensor") {
                                        model <- lm(y~P-1,weights=weights,singular.ok=singular.ok)
                                        coef.model <- coef(model)
                                        deriv.ind.vec[1:dim.P.tensor] <- TRUE
                                    } else if(basis.vec[p]=="taylor") {
                                        model <- lm(y~P,weights=weights,singular.ok=singular.ok)
                                        coef.model <- coef(model)[-1]
                                        deriv.ind.vec[1:dim.P.tensor] <- TRUE
                                    }
                                    model.deriv <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%coef.model[deriv.ind.vec]
                                } else {
                                    model.deriv <- rep(0, num.eval.obs)
                                }
                                
                            } else {
    
                                ## Compute factor `derivatives'
                                ## (differences), require base levels
    
                                P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,knots="quantiles",basis=basis.vec[p]))
                                if(attr(P,"relevant")) {
                                    if(basis.vec[p]=="additive" | basis.vec[p]=="taylor") {
                                        model.ma <- lm(y~P,weights=weights,singular.ok=singular.ok)
                                    } else {
                                        model.ma <- lm(y~P-1,weights=weights,singular.ok=singular.ok)
                                    }
                                    
                                    zeval.base.tmp <- zeval
                                    zeval.base.tmp[,kk] <- zeval.base[1,kk]
                                    
                                    P <- suppressWarnings(prod.spline(x=x,z=z,K=DS,I=include.vec,xeval=xeval,zeval=zeval.base.tmp,knots="quantiles",basis=basis.vec[p]))
                                    
                                    fit.spline <- suppressWarnings(predict(model.ma,newdata=data.frame(as.matrix(P))))
                                } else {
                                    model.ma <- lm(y~1,weights=weights,singular.ok=singular.ok)
                                    fit.spline <- suppressWarnings(predict(model.ma,newdata=data.frame(1:num.eval.obs)))
                                }
                                model.deriv <- fitted.mat[,p] - fit.spline
                            }
    
                        }
                        
                        deriv.mat[,p,k] <- model.deriv
    
                    }

                }
    
                list(fitted.mat=fitted.mat[,p],
                     rank.vec=rank.vec[p],
                     deriv.mat=deriv.mat[,p,])
                
            }
            
            for(p in 1:P.num) {
                fitted.mat[,p] <- output[[p]]$fitted.mat
                rank.vec[p] <- output[[p]]$rank.vec
                if(compute.deriv) deriv.mat[,p,] <- output[[p]]$deriv.mat
            }
            
            stopCluster(cl)
                    
        } ## end parallel
    
    }

    ## Compute fitted values and derivatives if requested
    
    if(verbose) cat("\r                                                                                               ")
    if(verbose) cat("\r")
    if(verbose) cat("\rComputing fitted values...")

    if(compute.mean) {
        fitted.values <- fitted.mat%*%b
    }

    if(compute.deriv) {
        if(verbose) cat("\r                                                                                               ")
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
    
    if(any(basis.singular.vec[b[b>ma.weights.cutoff]]) & verbose) warning("non-zero weight candidate model basis is ill-conditioned - reduce degree.max")

    if(verbose) {
        cat("\r                                                                                               ")
        cat("\r")
    }

    return(list(DKL.mat=if(is.null(ma.weights)){DKL.mat}else{DKL.mat.orig},
                S=S,
                X=X,
                all.combinations=all.combinations,
                auto.basis=auto.basis,
                auto.reduce=auto.reduce,
                auto.reduce.invoked=auto.reduce.invoked,
                basis.singular.vec=basis.singular.vec,
                basis.vec=if(is.null(ma.weights)){basis.vec}else{basis.vec.orig},
                basis=basis,
                eps.lambda=eps.lambda,
                compute.deriv=compute.deriv,
                compute.mean=compute.mean,
                degree.by=degree.by,
                degree.max=degree.max,
                degree.min=degree.min,
                deriv.index=deriv.index,
                deriv.order=deriv.order,
                deriv=deriv,
                knots=knots,
                lambda.S=lambda.S,
                lambda.num.max=lambda.num.max,
                lambda.seq=lambda.seq,
                ma.model.rank=sum(rank.vec*abs(b)),
                ma.weights=if(is.null(ma.weights)){abs(b)}else{ma.weights.orig},
                ma.weights.cutoff=ma.weights.cutoff,
                max.num.candidate.models=max.num.candidate.models,
                method=method,
                num.obs=num.obs,
                num.x=num.x,
                num.z=num.z,
                numeric.logical=numeric.logical,
                parallel.cores=parallel.cores,
                parallel=parallel,
                rank.vec=if(is.null(ma.weights)){rank.vec}else{rank.vec.orig},
                restrict.sum.ma.weights=restrict.sum.ma.weights,
                segments.by=segments.by,
                segments.max=segments.max,
                segments.min=segments.min,
                singular.ok=singular.ok,
                trace=trace,
                vc=vc,
                verbose=verbose,
                xnames=xnames,
                y=y,
                znames=znames,
                fitted.values=fitted.values))

}
