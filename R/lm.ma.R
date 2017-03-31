lm.ma <- function(...) UseMethod("lm.ma")

lm.ma.default <- function(y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("auto","tensor","glp","additive"),
                          compute.deriv=FALSE,
                          deriv.order=1,
                          degree.min=1,
                          degree.max=NULL,
                          lambda=1e-02,
                          segments.max=3,
                          knots=FALSE,
                          S=10,
                          method=c("jma","mma"),
                          ma.weights=NULL,
                          basis.vec=NULL,
                          bootstrap.ci=FALSE,
                          B=199,
                          alpha=0.05,
                          weights=NULL,
                          vc=TRUE,
                          verbose=TRUE,
                          tol=1e-08,
                          ...) {

    basis <- match.arg(basis)
    method <- match.arg(method)

    if(getOption("crs.messages")) {
        options(crs.messages=FALSE)
        exists.crs.messages <- TRUE
    } else {
        exists.crs.messages <- FALSE
    }
    
    if(!is.null(degree.max)) if(degree.max < 2) stop("You must average over at least two models")
    if(is.null(y) | is.null(X)) stop("You must provide data for y and X")
    if(!is.null(X) & !is.null(X.eval) & NCOL(X)!=NCOL(X.eval)) stop("X and X.eval must contain the same number of predictors")
    if(compute.deriv & (degree.min < deriv.order)) stop("Minimum degree (degree.min) must be at least as large the order of the derivative required (deriv.order)")
    if(degree.min < 1) stop("Minimum degree (degree.min) must be at least one")

    ## First obtain weights, then in subsequent call computes fits and 
    ## derivatives

    Est <- lm.ma.Est(y=y,
                     X=X,
                     X.eval=X.eval,
                     basis=basis,
                     compute.deriv=FALSE,
                     deriv.order=deriv.order,
                     degree.min=degree.min,
                     degree.max=degree.max,
                     lambda=lambda,
                     segments.max=segments.max,
                     knots=knots,
                     S=S,
                     method=method,
                     ma.weights=ma.weights,
                     basis.vec=basis.vec,
                     weights=weights,
                     vc=vc,
                     verbose=verbose,
                     tol=tol,
                     ...)
    
    ## Save rank vector and matrix of degrees and segments (not computed in next 
    ## call since ma.weights is passed, but needed for summary)
    
    if(compute.deriv | !is.null(X.eval)) {
        K.rank <- Est$rank.vec
        DS <- Est$DS
        basis.vec <- Est$basis.vec

        Est <- lm.ma.Est(y=y,
                         X=X,
                         X.eval=X.eval,
                         basis=basis,
                         compute.deriv=compute.deriv,
                         deriv.order=deriv.order,
                         degree.min=degree.min,
                         degree.max=degree.max,
                         lambda=lambda,
                         segments.max=segments.max,
                         knots=knots,
                         S=S,
                         method=method,
                         ma.weights=Est$ma.weights,
                         basis.vec=Est$basis.vec,
                         weights=weights,
                         vc=vc,
                         verbose=verbose,
                         tol=tol,
                         ...)
        
        Est$rank.vec <- K.rank
        Est$DS <- DS
        Est$DS <- basis.vec

    }

    ## Bootstrap is requested using pre-computed ma.weights
    
    if(bootstrap.ci) {

        options(crs.messages=FALSE)
    
        if(is.null(X.eval)) {
            n <- NROW(X) 
        } else {
            n <- NROW(X.eval)
        }
    
        boot.mat <- matrix(NA,n,B)
        if(compute.deriv) boot.deriv.array <- array(NA,c(n,B,Est$num.x))
    
        for(b in 1:B) {
            if(verbose) cat(paste("\rBootstrap replication ",b," of ",B,sep=""))
            ii <- sample(1:n,replace=TRUE)
            out.boot <- lm.ma.Est(y=y[ii],
                                  X=X[ii,],
                                  X.eval=if(is.null(X.eval)) { X } else {X.eval},
                                  basis=basis,
                                  compute.deriv=compute.deriv,
                                  deriv.order=deriv.order,
                                  degree.min=degree.min,
                                  degree.max=degree.max,
                                  lambda=lambda,
                                  segments.max=segments.max,
                                  knots=knots,
                                  S=S,
                                  method=method,
                                  ma.weights=Est$ma.weights,
                                  basis.vec=Est$basis.vec,
                                  weights=weights,
                                  vc=vc,
                                  verbose=FALSE,
                                  tol=1e-08,
                                  ...)
            boot.mat[,b] <- out.boot$fitted
            if(compute.deriv) for(k in 1:Est$num.x) boot.deriv.array[,b,k] <- out.boot$deriv[,k]
        }

        if(verbose) cat("\r                                     ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing quantiles...")

        Est$fitted.ci.l <- apply(boot.mat,1,quantile,prob=alpha/2,type=1,na.rm=TRUE)
        Est$fitted.ci.u <- apply(boot.mat,1,quantile,prob=1-alpha/2,type=1,na.rm=TRUE)
        
        if(compute.deriv) {
            Est$deriv.ci.l <- Est$deriv.ci.u <- matrix(NA,n,Est$num.x)
            for(k in 1:Est$num.x) {
                Est$deriv.ci.l[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=alpha/2,type=1,na.rm=TRUE)
                Est$deriv.ci.u[,k] <- apply(boot.deriv.array[,,k],1,quantile,prob=1-alpha/2,type=1,na.rm=TRUE)        
            }
        }
        if(verbose) cat("\r                                     ")
        if(verbose) cat("\r")
    }

    if(exists.crs.messages) options(crs.messages=TRUE)

    class(Est) <- "lm.ma"
    return(Est)

}

lm.ma.Est <- function(y=NULL,
                      X=NULL,
                      X.eval=NULL,
                      basis=c("tensor","glp","additive","auto"),
                      compute.deriv=FALSE,
                      deriv.order=1,
                      degree.min=1,
                      degree.max=NULL,
                      lambda=1e-02,
                      segments.max=3,
                      knots=FALSE,
                      S=10,
                      method=c("jma","mma"),
                      ma.weights=NULL,
                      basis.vec=NULL,
                      weights=NULL,
                      vc=TRUE,
                      verbose=TRUE,
                      tol=1e-08,
                      ...) {
    
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
        xeval <- x
        zeval <- z
    }

    rm(xztmp)

    if(vc & !is.null(num.z)) {
        z.unique <- crs:::uniquecombs(as.matrix(z))
        num.z <- ncol(z.unique)
        ind <-  attr(z.unique,"index")
        ind.vals <-  unique(ind)
        nrow.z.unique <- nrow(z.unique)
        zeval.unique <- crs:::uniquecombs(as.matrix(zeval))
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
        degree.max <- max(3,floor((S/num.x)*(NROW(X)/100)^0.1))
    }
    
    P <- degree.max
    ma.weights.orig <- ma.weights
    basis.vec.orig <- basis.vec

    if(knots) {
        K.mat <- crs:::matrix.combn(K.vec1=degree.min:degree.max,K.vec2=1:segments.max,num.x=num.x)
    } else {
        K.mat <- crs:::matrix.combn(K.vec1=degree.min:degree.max,num.x=num.x)
        K.mat <- cbind(K.mat[,1:num.x],matrix(1,nrow(K.mat),num.x,byrow=TRUE))
    }
    if(basis=="auto" & is.null(basis.vec) & is.null(ma.weights)) {
        basis.vec <- character()
    } else if(basis=="auto" & !is.null(basis.vec) & !is.null(ma.weights)) {
        basis.vec <- basis.vec[ma.weights>1e-05]
    }  else if(basis!="auto") {
        basis.vec <- rep(basis,nrow(K.mat))
    }
        if(!is.null(ma.weights)) {
        K.mat <- K.mat[ma.weights>1e-05,]
        ma.weights <- ma.weights[ma.weights>1e-05]/sum(ma.weights[ma.weights>1e-05])
    }

    P.num <- NROW(K.mat)

    deriv <- NULL
    if(compute.deriv) {
        deriv.mat <- array(NA,c(if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},P.num,num.x))
        deriv <- matrix(NA,if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},num.x)
        colnames(deriv) <- xnames
    }

    K.rank <- numeric(length=P.num)
    sigsq <- numeric(length=P.num)
    ma.mat <- matrix(NA,NROW(X),P.num)
    fitted.mat <- matrix(NA,if(is.null(X.eval)){NROW(X)}else{NROW(X.eval)},P.num)

    for(p in P.num:1) {
        
        DS <- cbind(K.mat[p,1:num.x],K.mat[p,(num.x+1):(2*num.x)])   

        if(verbose) cat(paste("\rCandidate model ",P.num-p+1," of ",P.num," (degree.max = ",degree.max,")...",sep=""))

        if(is.null(ma.weights)) {

            if(vc & !is.null(num.z)) {

                if(basis=="auto") {
                    cv.min <- Inf
                    for(b.basis in c("tensor","glp","additive")) {

                        fit.spline <- numeric(length=NROW(x))
                        htt <- numeric(length=NROW(x))
                        basis.singular <- logical(length=nrow.z.unique)
                        for(i in 1:nrow.z.unique) {
                            zz <- ind == ind.vals[i]
                            L <- crs:::prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                            if(!is.null(weights)) L <- weights*L
                            P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=b.basis)
                            if(!is.fullrank(P)) basis.singular[i] <- TRUE
                            if(b.basis=="additive" || b.basis=="glp") {
                                model.z.unique <- lm(y~P,weights=L)
                            } else {
                                model.z.unique <- lm(y~P-1,weights=L)
                            }
                            htt[zz] <- hatvalues(model.z.unique)[zz]
                            P <- suppressWarnings(crs:::prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=b.basis))
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

                    if(basis.singular.min) warning("dimension basis is ill-conditioned - reduce degree.max")

                    fit.spline <- fit.spline.min
                    htt.min <- htt
                    model.z.unique$rank <- rank.min

                } else {

                    fit.spline <- numeric(length=NROW(x))
                    htt <- numeric(length=NROW(x))
                    basis.singular <- logical(length=nrow.z.unique)
                    for(i in 1:nrow.z.unique) {
                        zz <- ind == ind.vals[i]
                        L <- crs:::prod.kernel(Z=z,z=z.unique[ind.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z)
                        if(!is.null(weights)) L <- weights*L
                        P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                        if(!is.fullrank(P)) basis.singular[i] <- TRUE
                        if(basis.vec[p]=="additive" || basis.vec[p]=="glp") {
                            model.z.unique <- lm(y~P,weights=L)
                        } else {
                            model.z.unique <- lm(y~P-1,weights=L)
                        }
                        htt[zz] <- hatvalues(model.z.unique)[zz]
                        P <- suppressWarnings(crs:::prod.spline(x=x,K=DS,xeval=x[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                        fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))
                    }
                    if(any(basis.singular==TRUE)) warning("dimension basis is ill-conditioned - reduce degree.max")
                 }

                fitted.mat[,p] <- fit.spline
                    
                if(method=="mma") {
                    ma.mat[,p] <- y - fit.spline
                } else {
                    htt <- ifelse(htt < 1, htt, 1-.Machine$double.eps)
                    ma.mat[,p] <- fit.spline - htt*(y - fit.spline)/(1-htt)
                }

                K.rank[p] <- model.z.unique$rank
                sigsq[p] <- sqrt(sum((y - fit.spline)^2)/(NROW(x)-model.z.unique$rank))

            } else {

                if(basis=="auto") {
                    cv.min <- Inf
                    for(b.basis in c("tensor","glp","additive")) {
                        P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=b.basis)
                        if(!is.fullrank(P)) basis.singular <- TRUE
                        if(b.basis=="additive" || b.basis=="glp") {
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
                    
                    if(basis.singular.min) warning("dimension basis is ill-conditioned - reduce degree.max")
                    
                    model.ma <- model.ma.min
                    fit.spline <- fitted(model.ma)

                } else {

                    P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                    if(!is.fullrank(P)) warning("dimension basis is ill-conditioned - reduce degree.max")
                    if(basis.vec[p]=="additive" || basis.vec[p]=="glp") {
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

                K.rank[p] <- model.ma$rank

                sigsq[p] <- sqrt(sum(residuals(model.ma)^2)/(NROW(x)-model.ma$rank))

            }
            
        }

        if(!is.null(ma.weights))  {
            if(vc & !is.null(num.z)) {

                fit.spline <- numeric(length=num.eval)
                for(i in 1:nrow.zeval.unique) {
                    zz <- ind.zeval == ind.zeval.vals[i]
                    L <- suppressWarnings(crs:::prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z))
                    if(!is.null(weights)) L <- weights*L
                    P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                    if(basis.vec[p]=="additive" || basis.vec[p]=="glp") {
                        model.z.unique <- lm(y~P,weights=L)
                    } else {
                        model.z.unique <- lm(y~P-1,weights=L)
                    }
                    P <- suppressWarnings(crs:::prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p]))
                    fit.spline[zz] <- suppressWarnings(predict(model.z.unique,newdata=data.frame(as.matrix(P))))

                }

                fitted.mat[,p] <- fit.spline

            } else {

                P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                if(basis.vec[p]=="additive" || basis.vec[p]=="glp") {
                    model.ma <- lm(y~P,weights=weights)
                } else {
                    model.ma <- lm(y~P-1,weights=weights)
                }
                
                P <- suppressWarnings(crs:::prod.spline(x=x,z=z,K=DS,I=include,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p]))

                fitted.mat[,p] <- suppressWarnings(predict(model.ma,newdata=data.frame(as.matrix(P))))
                
            }
            
            if(compute.deriv) {

                if(basis.vec[p]=="additive" || basis.vec[p]=="glp") {
                    K.additive <- DS
                    K.additive[,2] <- ifelse(DS[,1]==0,0,DS[,2])
                    K.additive[,1] <- ifelse(DS[,1]>0,DS[,1]-1,DS[,1])
                }
                
                for(k in 1:num.x) {
                   
                    if(vc & !is.null(num.z)) {

                        deriv.spline <- numeric(length(num.eval))
                        for(i in 1:nrow.zeval.unique) {
                            zz <- ind.zeval == ind.zeval.vals[i]
                            L <- suppressWarnings(crs:::prod.kernel(Z=z,z=zeval.unique[ind.zeval.vals[i],],lambda=lambda.vec,is.ordered.z=is.ordered.z))
                            if(!is.null(weights)) L <- weights*L
                            P <- crs:::prod.spline(x=x,K=DS,knots="quantiles",basis=basis.vec[p])
                            P.deriv <- suppressWarnings(crs:::prod.spline(x=x,K=DS,xeval=xeval[zz,,drop=FALSE],knots="quantiles",basis=basis.vec[p],deriv.index=k,deriv=deriv.order))
                            if(basis.vec[p]=="additive") {
                                model <- lm(y~P,weights=L)
                                dim.P.deriv <- sum(K.additive[k,])
                                deriv.start <- ifelse(k!=1,sum(K.additive[1:(k-1),])+1,1)
                                deriv.end <- deriv.start+sum(K.additive[k,])-1
                                deriv.ind.vec <- deriv.start:deriv.end
                                deriv.spline[zz] <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%(coef(model)[-1])[deriv.ind.vec]
                            } else if(basis.vec[p]=="tensor") {
                                model <- lm(y~P-1,weights=L)
                                deriv.spline[zz] <- P.deriv%*%coef(model)
                            } else if(basis.vec[p]=="glp") {
                                model <- lm(y~P,weights=L)
                                deriv.spline[zz] <- P.deriv%*%coef(model)[-1]
                            }
                        }

                        model.deriv <- deriv.spline

                    } else {

                        P <- crs:::prod.spline(x=x,z=z,K=DS,I=include,knots="quantiles",basis=basis.vec[p])
                        P.deriv <- suppressWarnings(crs:::prod.spline(x=x,z=z,K=DS,I=include,xeval=xeval,zeval=zeval,knots="quantiles",basis=basis.vec[p],deriv.index=k,deriv=deriv.order))
                        dim.P.tensor <- NCOL(P)
                        deriv.ind.vec <- logical(length=NCOL(P))
                        coef.vec.model <- numeric(length=NCOL(P))
                        if(basis.vec[p]=="additive") {
                            model <- lm(y~P,weights=weights)
                            coef.vec.model <- coef(model)[-1]
                            dim.P.deriv <- sum(K.additive[k,])
                            deriv.start <- ifelse(k!=1,sum(K.additive[1:(k-1),])+1,1)
                            deriv.end <- deriv.start+sum(K.additive[k,])-1
                            deriv.ind.vec[deriv.start:deriv.end] <- TRUE
                        } else if(basis.vec[p]=="tensor") {
                            model <- lm(y~P-1,weights=weights)
                            coef.vec.model <- coef(model)
                            deriv.ind.vec[1:dim.P.tensor] <- TRUE
                        } else if(basis.vec[p]=="glp") {
                            model <- lm(y~P,weights=weights)
                            coef.vec.model <- coef(model)[-1]
                            deriv.ind.vec[1:dim.P.tensor] <- TRUE
                        }
                        
                        model.deriv <- P.deriv[,deriv.ind.vec,drop=FALSE]%*%coef.vec.model[deriv.ind.vec]

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

        ## Solve the quadratic program for the Mallows model average weights
        M <- ncol(ma.mat)
        D <- t(ma.mat)%*%ma.mat
        #if(qr(D)$rank<M) D <- D + diag(tol,M,M)
        tol.ridge <- tol
        #while(!is.fullrank(D)) {
        while(qr(D)$rank<M) {
            D <- D + diag(tol.ridge,M,M)
            tol.ridge <- tol.ridge*10
            print("ridging")
            print(tol.ridge)
        }
        A <- cbind(rep(1,M),diag(1,M,M))
        b0 <- c(1,rep(0,M))
        if(method=="mma") {
            d <- -sigsq[which.max(K.rank)]*K.rank
        } else {
            d <- t(y)%*%ma.mat
        }        
        b <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution

        rebalance <- TRUE
        if(rebalance) {
            ## Solve the quadratic program for the Mallows model average weights
            ma.mat.reb <- ma.mat[,b>1e-05,drop=FALSE]
            M <- ncol(ma.mat.reb)
            D <- t(ma.mat.reb)%*%ma.mat.reb
            #if(qr(D)$rank<M) D <- D + diag(tol,M,M)
            tol.ridge <- tol
            #while(!is.fullrank(D)) {
            while(qr(D)$rank<M) {
                D <- D + diag(tol.ridge,M,M)
                tol.ridge <- tol.ridge*10
                print("rebalance ridging")
                print(tol.ridge)
            }
            A <- cbind(rep(1,M),diag(1,M,M))
            b0 <- c(1,rep(0,M))
            K.rank.reb <- K.rank[b>1e-05]
            if(method=="mma") {
                d <- -sigsq[which.max(K.rank)]*K.rank.reb
            } else {
                d <- t(y)%*%ma.mat.reb
            }        
            b.reb <- solve.QP(Dmat=D,dvec=d,Amat=A,bvec=b0,meq=1)$solution

            if(!isTRUE(all.equal(as.numeric(b[b>1e-05]),as.numeric(b.reb)))) {
                warning("rebalancing occured")
                b[b>1e-05] <- b.reb
            }
        }

    } else {
        ## For bootstrapping use weights from initial call
        b <- ma.weights
    }

    if(compute.deriv) {
        if(verbose) cat("\r                                                    ")
        if(verbose) cat("\r")
        if(verbose) cat("\rComputing derivatives...")
        for(k in 1:num.x) deriv[,k] <- deriv.mat[,,k]%*%b
    }

    if(verbose) {
        cat("\r                                                    ")
        cat("\r")
    }

    return(list(fitted.values=fitted.mat%*%b,
                deriv=deriv,
                ma.weights=if(is.null(ma.weights)){abs(b)}else{ma.weights.orig},
                basis.vec=if(is.null(ma.weights)){basis.vec}else{basis.vec.orig},
                y=y,
                X=X,
                basis=basis,
                compute.deriv=compute.deriv,
                deriv.order=deriv.order,
                degree.min=degree.min,
                degree.max=degree.max,
                lambda=lambda,
                segments.max=segments.max,
                knots=knots,
                S=S,
                method=method,
                num.x=num.x,
                num.z=num.z,
                rank.vec=K.rank,
                nobs=NROW(y),
                DS=K.mat,
                vc=vc,
                verbose=verbose,
                tol=tol))

}

lm.ma.formula <- function(formula,
                          data=list(),
                          y=NULL,
                          X=NULL,
                          X.eval=NULL,
                          basis=c("auto","tensor","glp","additive"),
                          compute.deriv=FALSE,
                          deriv.order=1,
                          degree.min=1,
                          degree.max=NULL,
                          lambda=1e-02,
                          segments.max=3,
                          knots=FALSE,
                          S=10,
                          method=c("jma","mma"),
                          ma.weights=NULL,
                          basis.vec=NULL,
                          bootstrap.ci=FALSE,
                          B=199,
                          alpha=0.05,
                          weights=NULL,
                          vc=TRUE,
                          verbose=TRUE,
                          tol=1e-08,
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
                       degree.max=degree.max,
                       lambda=lambda,
                       segments.max=segments.max,
                       knots=knots,
                       S=S,
                       method=method,
                       ma.weights=ma.weights,
                       basis.vec=basis.vec,
                       bootstrap.ci=bootstrap.ci,
                       B=B,
                       alpha=alpha,
                       weights=weights,
                       vc=vc,
                       verbose=verbose,
                       tol=tol,
                       ...)

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
  cat("\nModel Averaging Linear Regression",sep="")
  cat(paste(ifelse(object$vc, " (Varying Coefficient Specification)"," (Additive Dummy Specification)"),sep=""))
  cat(paste("\nModel average criterion: ", ifelse(object$method=="jma","Jacknife (Hansen and Racine (2013))","Mallows  (Hansen (2007))"), sep=""))
  cat(paste("\nMinimum degree: ", object$degree.min, sep=""))  
  cat(paste("\nMaximum degree: ", object$degree.max, sep=""))
  cat(paste("\nBasis: ", object$basis, sep=""))  
  cat(paste("\nNumber of observations: ", object$nobs, sep=""))
  cat(paste("\nEquivalent number of parameters: ", formatC(sum(object$rank.vec*object$ma.weights),format="f",digits=2), sep=""))
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
                         lambda=object$lambda,
                         segments.max=object$segments.max,
                         knots=object$knots,
                         S=object$S,
                         method=object$method,
                         ma.weights=object$ma.weights,
                         basis.vec=object$basis.vec,
                         weights=object$weights,
                         vc=object$vc,
                         verbose=object$verbose,
                         tol=object$tol,
                         ...)
    fitted.values <- Est$fitted.values
    deriv <- Est$deriv
  }

  if(is.null(Est$deriv)) {
      return(fitted.values)
  } else {
      return(list(fit=fitted.values,deriv=deriv))
  }

}
