
#############################################################
## Type 1 tensor product methods start here (i.e. Wood, 2006)
#############################################################

tensor.prod.model.matrix <- function(X) {
# X is a list of model matrices, from which a tensor product model matrix is to be produced.
# e.g. ith row is basically X[[1]][i,]%x%X[[2]][i,]%x%X[[3]][i,], but this routine works 
# column-wise, for efficiency, and does work in compiled code.
  m <- length(X)              ## number to row tensor product
  d <- unlist(lapply(X,ncol)) ## dimensions of each X
  n <- nrow(X[[1]])           ## columns in each X
  X <- as.numeric(unlist(X))  ## append X[[i]]s columnwise
  T <- numeric(n*prod(d))     ## storage for result
  .C("mgcv_tensor_mm",
     as.double(X),
     as.double(T),
     as.integer(d),
     as.integer(m),
     as.integer(n),
     package="ma") ## produce product
  ## Give T attributes of matrix. Note that initializing T as a matrix 
  ## requires more time than forming the row tensor product itself (R 3.0.1)
  attr(T,"dim") <- c(n,prod(d)) 
  class(T) <- "matrix"
  T
} ## end tensor.prod.model.matrix

uniquecombs<-function(x) {
## takes matrix x and counts up unique rows
## `unique' now does this in R
if (is.null(x)) stop("x is null")
if (is.null(nrow(x))) stop("x has no row attribute")
if (is.null(ncol(x))) stop("x has no col attribute")
ind <- rep(0,nrow(x))
res<-.C("RuniqueCombs",x=as.double(x),ind=as.integer(ind),
        r=as.integer(nrow(x)),c=as.integer(ncol(x)))
n <- res$r*res$c
x <- matrix(res$x[1:n],res$r,res$c)
attr(x,"index") <- res$ind+1 ## C to R index gotcha
x
}

