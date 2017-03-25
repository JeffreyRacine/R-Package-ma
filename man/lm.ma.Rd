\name{lm.ma}
\alias{lm.ma}
\alias{lm.ma.default}
\alias{lm.ma.formula}
\title{
Model Averaging
}
\description{
A function similar to \code{\link{lm}} that averages over a set of candidate models.
}
\usage{
lm.ma(...)

\method{lm.ma}{default}(y=NULL,
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
         ...)

\method{lm.ma}{formula}(formula,
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
         ...)

}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit }

  \item{data}{ an optional data frame containing the variables in the model }

  \item{y}{a one dimensional vector of dependent data}

  \item{X}{a \eqn{p}-variate data frame of explanatory (training) data }

  \item{X.eval}{a \eqn{p}-variate data frame of points on which the regression will be estimated (evaluation data)}

\item{basis}{ a character string indicating whether the additive or tensor product basis matrix for a multivariate polynomial should be used}

  \item{compute.deriv}{a logical value indicating whether to compute derivatives or not}

  \item{p.max}{maximum value for the basis in each dimension (defaults to \code{round((10/NCOL(X))*(NROW(X)/100)^0.25)}}

  \item{method}{whether to use jackknife model averaging (JMA) or Mallows model averaging (MMA)}

  \item{ma.weights}{a vector of weights obtained from a previous invocation (useful for bootstrapping)}
  
  \item{bootstrap.ci}{a logical value indicating whether to bootstrap nonparametric confidence intervals or not}
  
  \item{B}{the number of bootstrap replications desired}
  
  \item{alpha}{value in (0,1) used to compute 1-alpha\% confidence intervals}
  
  \item{...}{optional arguments to be passed}

}
\details{
Yes...
}
\value{

  \item{fitted.values}{vector of fitted values}
  \item{deriv}{matrix of derivative vectors for each predictor}
  \item{ma.weights}{model average weights}
  \item{p.max}{maximum value of p.max for each dimension}
  \item{fitted.ci.l}{alpha/2 nonparametric confidence value for the fitted/predicted values}
  \item{fitted.ci.u}{1-alpha/2 nonparametric confidence value for the fitted/predicted values}  
  \item{deriv.ci.l}{alpha/2 nonparametric confidence value matrix for
    the derivatives}
  \item{deriv.ci.u}{1-alpha/2 nonparametric confidence value matrix for
    the derivatives}
  \item{r.squared}{nonlinear measure of goodness of fit}
  \item{residuals}{model residuals}

}
\references{
YO
}
\author{
Jeffrey S. Racine
}
\note{
This codes is in beta status until further notice - proceed accordingly.
}

\seealso{
\code{\link{lm}}
}
\examples{
require(rgl)

set.seed(42)
n <- 1000
x1 <- runif(n,-1,1)
x2 <- runif(n,-1,1)

#dgp <- x1+x2
#dgp <- x1**2+x2**2

dgp <- cos(pi*x1)*sin(pi*x2)

y <- dgp + rnorm(n,sd=0.5)

n.eval <- 25
x.seq <- seq(-1,1,length=n.eval)
newdata <- data.frame(expand.grid(x.seq,x.seq))
names(newdata) <- c("x1","x2")

model <- lm.ma(y~x1+x2)

z <- matrix(predict(model,newdata=newdata),n.eval,n.eval)

num.colors <- 1000
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

## Open an rgl 3d window and use `persp3d()', a high-level function
## for 3D surfaces (and define the size of the window to be
## 640x640). The function par3d() passes in a window size (the default
## is 256x256 which is quite small), the function rgl.viewpoint()
## allows you to modify the `field of view' to get more of a
## `perspective' feel to the plot, while the function grid3d() adds a
## grid to the plot.

open3d()

par3d(windowRect=c(900,100,900+640,100+640))
rgl.viewpoint(theta = 0, phi = -70, fov = 80)

persp3d(x.seq,x.seq,z=z,
        xlab="X1",ylab="X2",zlab="Y",
        ticktype="detailed",
        border="red",
        color=col,
        alpha=.7,
        back="lines",
        main="Conditional Mean")

grid3d(c("x", "y+", "z"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
