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
    quantile(x, probs = prob)
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
