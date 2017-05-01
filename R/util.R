## This function provides a matrix with all combinations of a vector K
## (containing degrees) and vector lambda.

## Always must have knots and segments

matrix.combn <- function(K.vec1,K.vec2=NULL,K.vec3=NULL,num.x=0,num.z=0) {

  if(num.x==0 & num.z==0) stop(" must provide at least one variable")

  ls <- list()
  for(i in 1:num.x) ls[[i]] <- K.vec1
  if(!is.null(K.vec2)) for(i in 1:num.x) ls[[num.x+i]] <- K.vec2
  if(!is.null(K.vec3)) for(i in 1:num.z) ls[[2*num.x+i]] <- K.vec3
  return(as.matrix(do.call(expand.grid,ls)))

}

lambda.plugin <- function(z) {

    ## Optimal bandwidth for univariate probability function estimation...

    n <- length(z)
    p.z <- prop.table(table(z))
    C <- length(p.z)
    Lambda.1 <- (1-C*p.z)/(C-1)
    Lambda.2 <- (1+C**2*p.z - 2*C*p.z)/(C-1)**2
    Lambda.3 <- p.z*(1+Lambda.1)
    return(sum(Lambda.2)/(sum(Lambda.2-Lambda.1**2) + n*sum(Lambda.1**2)))

}

blank <- function(len){
  sapply(len, function(nb){
    paste(rep(' ', times = nb), collapse='')
  })
}

uocquantile = function(x, prob = 0.5) {
  if (is.ordered(x)){
    tq = unclass(table(x))
    tq = tq / sum(tq)
    j = which(sapply(1:length(tq), function(y){ sum(tq[1:y]) }) >= prob)[1]
    sort(unique(x))[j]
  } else if (is.factor(x)) {
    ## just returns mode
    tq = unclass(table(x))
    j = which(tq == max(tq))[1]
    sort(unique(x))[j]
  } else {
    quantile(x, probs = prob, type = 1)
  }
}

is.fullrank <- function(x)
{
  e <- eigen(crossprod(as.matrix(x)), symmetric = TRUE, only.values = TRUE)$values
  e[1] > 0 && abs(e[length(e)]/e[1]) > max(dim(x))*max(sqrt(abs(e)))*.Machine$double.eps
}

splitFrame <- function(xz, factor.to.numeric=FALSE) {
  
  if(missing(xz)) stop(" you must provide xz data")
  if(!is.data.frame(xz)) stop(" xz must be a data frame")

  xznames <- names(xz)
  
  IND <- logical()

  for(i in 1:NCOL(xz)) IND[i] <- is.factor(xz[,i])

  x <- xz[,!IND,drop=FALSE]
  num.x <- ncol(x)

  ## We require at least one continuous predictor to conduct spline
  ## smoothing, but there may/may not be factors.

  if(num.x == 0) stop(" can't fit spline surfaces with no continuous predictors")

  xnames <- xznames[!IND]

  is.ordered.z <- NULL
  
  if(any(IND)) {
    is.ordered.z <- logical()
    for(i in 1:NCOL(xz[,IND,drop=FALSE])) is.ordered.z[i] <- is.ordered((xz[,IND,drop=FALSE])[,i])
    if(!factor.to.numeric) {
      z <- data.frame(xz[,IND,drop=FALSE])
    } else {
      ## If factor.to.numeric crudely convert factors to numeric.
      z <- matrix(NA,NROW(xz),NCOL(xz[,IND,drop=FALSE]))
      ## To "revert" a factor f to its original numeric values,
      ## as.numeric(levels(f))[f] is recommended. Problem is that for
      ## character strings it produces a warning message. No idea how
      ## to test for this so dropping for the moment. Will affect
      ## ordered types.
      for(i in 1:NCOL(xz[,IND,drop=FALSE])) {
        suppressWarnings(z[,i] <- as.numeric(levels((xz[,IND,drop=FALSE])[,i]))[(xz[,IND,drop=FALSE])[,i]])
        if(any(is.na(z[,i]))) z[,i] <- as.numeric((xz[,IND,drop=FALSE])[,i])
      }
    }
    ## Don't assign names when factor.to.numeric is TRUE (otherwise
    ## NAs populate matrix)
    if(!factor.to.numeric) names(z) <- xznames[IND]
    znames <- xznames[IND]
    num.z <- ncol(z)
  } else {
    z <- NULL
    znames <- NULL
    num.z <- NULL
  }

  return(list(x=x,
              num.x=num.x,
              xnames=xnames,
              z=z,
              num.z=num.z,
              numeric.logical=!IND,
              is.ordered.z=is.ordered.z,
              znames=znames))
  
}

RSQfunc <- function(y,y.pred,weights=NULL) {
  if(!is.null(weights)) {
    y <- y*sqrt(weights)
    y.pred <- y.pred*sqrt(weights)
  }
  y.mean <- mean(y)
  return((sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2)))
}

dim.bs <- function(basis="additive",kernel=TRUE,degree=NULL,segments=NULL,include=NULL,categories=NULL) {

  ## This function computes the dimension of the taylor basis without the
  ## memory overhead associated with computing the taylor basis itself
  ## (thanks to Zhenghua Nie)

  two.dimen<- function(d1,d2,nd1,pd12){
    if(d2 ==1) {
      ret <- list()
      ret$d12 <- pd12
      ret$nd1 <- nd1
      return(ret)
    }
    d12 <- d2
    if(d1-d2>0){
      for(i in 1:(d1-d2)){
        d12 <- d12+d2*nd1[i]
      }}
    if(d2>1){
      for(i in 2:d2){
        d12 <- d12 + (i*nd1[d1-i+1])
      }
    }
    d12 <- d12 + nd1[d1]   ## The maximum number
    
    nd2 <- nd1  ## Calculate nd2
    if(d1>1){
      for(j in 1:(d1-1)) {
        nd2[j] <- 0
        for(i in j:max(0,j-d2+1)) {
          if(i > 0) {
            nd2[j] <- nd2[j] + nd1[i]                  
          }
          else {
            nd2[j] <- nd2[j] + 1  ## nd1[0] always 1
          }
        }
      }
    }
    if(d2>1) {
      nd2[d1] <- nd1[d1]
      for(i in (d1-d2+1):(d1-1)) nd2[d1] <- nd2[d1]+nd1[i]
    }
    else {
      nd2[d1] <- nd1[d1]
    }
    ret <- list()
    ret$d12 <- d12
    ret$nd1 <- nd2 
    
    return(ret)
  }
  
  ## Some basic error checking
 
  if(basis!="additive" & basis!="taylor" & basis!="tensor") stop(" Error: basis must be either additive, taylor, or tensor")

  #if(!kernel)
  #  if(is.null(include) | is.null(categories)) stop(" Error: you must provide include and categories vectors")    
  
  K <- cbind(degree,segments)

  ncol.bs <- 0

  if(kernel) {
    if(basis=="additive") {
      if(any(K[,1] > 0))
        ncol.bs <- sum(rowSums(K[K[,1]!=0,,drop=FALSE])-1)
    }
    if(basis=="taylor") {
      dimen <- rowSums(K[K[,1]!=0,,drop=FALSE])-1
      dimen <- dimen[dimen>0] ## Delete elements which are equal to 0.
      dimen <- sort(dimen,decreasing=TRUE) ## Sort the array to save memory when doing the computation.
      k <-length(dimen)
      if(k==0) {
        ncol.bs <- 0
      } else {
        nd1 <- rep(1,dimen[1])   ## At the beginning,  we have one for [1, 2, 3, ..., dimen[1]]
        nd1[dimen[1]] <- 0       ## nd1 represents the frequency for every element of [1, 2, 3, ..., dimen[1]]
        ncol.bs <- dimen[1]
        if(k>1) {
          for(i in 2:k) {
            dim.rt <- two.dimen(dimen[1],dimen[i],nd1,ncol.bs)
            nd1 <- dim.rt$nd1
            ncol.bs <- dim.rt$d12
          }
          ncol.bs <- dim.rt$d12+k-1
        }
      }
    }
    if(basis=="tensor") {
      if(any(K[,1] > 0))
        ncol.bs <- prod(rowSums(K[K[,1]!=0,,drop=FALSE]))
    }
  } else {
    if(basis=="additive") {
      if(any(K[,1] > 0)) 
        ncol.bs <- sum(c(rowSums(K[K[,1]!=0,,drop=FALSE]),include*categories-1))
    }
    if(basis=="taylor") {
      dimen <- c(rowSums(K[K[,1]!=0,,drop=FALSE])-1,include*categories-1)
      dimen <- dimen[dimen>0] ## Delete elements which are eqaul to 0.
      dimen <- sort(dimen,decreasing=TRUE) ## Sort the array to save memory when doing the computation.
      k <-length(dimen)
      if(k==0) {
        ncol.bs <- 0
      } else {
        nd1 <- rep(1,dimen[1])   ## At the beginning,  we have one for [1, 2, 3, ..., dimen[1]]
        nd1[dimen[1]] <- 0       ## nd1 represents the frequency for every element of [1, 2, 3, ..., dimen[1]]
        ncol.bs <- dimen[1]
        if(k>1) {
          for(i in 2:k) {
            dim.rt <- two.dimen(dimen[1],dimen[i],nd1,ncol.bs)
            nd1 <- dim.rt$nd1
            ncol.bs <- dim.rt$d12
          }
          ncol.bs <- dim.rt$d12+k-1
        }
      }
    }
    if(basis=="tensor") {
      if(any(K[,1] > 0)) 
        ncol.bs <- prod(c(rowSums(K[K[,1]!=0,,drop=FALSE]),(include*categories-1)))
    }
  }

  return(ncol.bs)

}

