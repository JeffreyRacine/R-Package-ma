\name{lm.ma}
\alias{lm.ma}
\alias{lm.ma.default}
\alias{lm.ma.formula}
\title{
Fitting Linear Model Average Models
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
      deriv.order=1,
      p.max=NULL,
      method=c("jma","mma"),
      ma.weights=NULL,
      bootstrap.ci=FALSE,
      B=100,
      alpha=0.05,
      weights=NULL,
      ...)

\method{lm.ma}{formula}(formula,
      data=list(),
      y=NULL,
      X=NULL,
      X.eval=NULL,
      basis=c("glp","tensor","additive"),
      compute.deriv=FALSE,
      deriv.order=1,
      p.max=NULL,
      method=c("jma","mma"),
      ma.weights=NULL,
      bootstrap.ci=FALSE,
      B=100,
      alpha=0.05,
      weights=NULL,
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
  
  \item{deriv.order}{integer indicating order of derivative desired}

  \item{p.max}{maximum value for the basis in each dimension (defaults to \code{round((10/NCOL(X))*(NROW(X)/100)^0.25)}}

  \item{method}{whether to use jackknife model averaging (jma) or Mallows model averaging (mma)}

  \item{ma.weights}{a vector of weights obtained from a previous invocation (useful for bootstrapping)}
  
  \item{bootstrap.ci}{a logical value indicating whether to bootstrap nonparametric confidence intervals or not}
  
  \item{B}{the number of bootstrap replications desired}
  
  \item{alpha}{value in (0,1) used to compute 1-alpha\% confidence intervals}
  
  \item{weights}{an optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector; if non-NULL, weighted least squares is used with weights \code{weights} (that is, minimizing \eqn{sum(w*e^2)}); otherwise ordinary least squares is used}
  
  \item{...}{optional arguments to be passed}

}
\details{
Models for lm.ma are specified symbolically. A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response.
}
\value{

lm returns an object of class "lm.ma".

The function summary is used to obtain and print a summary of the results. The generic accessor functions fitted.values and residuals extract various useful features of the value returned by lm.ma.

An object of class "lm.ma" is a list containing at least the following components:

  \item{fitted.values}{vector of fitted values}
  \item{deriv}{matrix of derivative vectors for each predictor}
  \item{ma.weights}{model average weights}
  \item{p.max}{value of p.max for each dimension}
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
TBD
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
## Example 1

set.seed(42)
n <- 100
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + rnorm(n,sd=0.5*sd(dgp))

model <- lm.ma(y~x,bootstrap.ci=TRUE,compute.deriv=TRUE)
model$ma.weights

par(mfrow=c(1,2))

plot(x,y,cex=0.25)
lines(x,fitted(model))
lines(x,model$fitted.ci.l,col=2,lty=2)
lines(x,model$fitted.ci.u,col=2,lty=2)

ylim <- range(c(model$deriv.ci.l[,1],model$deriv.ci.u[,1]))
plot(x,model$deriv[,1],ylim=ylim,ylab="derivative",type="l")
lines(x,model$deriv.ci.l[,1],col=2,lty=2)
lines(x,model$deriv.ci.u[,1],col=2,lty=2)

par(mfrow=c(1,1))

## Example 2

require(rgl)

set.seed(42)
n <- 1000
x1 <- runif(n)
x2 <- runif(n)

dgp <- cos(2*pi*x1)*sin(2*pi*x2)

y <- dgp + rnorm(n,sd=0.5*sd(dgp))

n.eval <- 25
x.seq <- seq(0,1,length=n.eval)
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

## Example 3

require(crs)
data(wage1)

## Classical linear regression model

model.lm <- lm(lwage ~ factor(female) + factor(married) + educ + exper + tenure, data = wage1)
summary(model.lm)$r.squared

## Compare with the model average estimator (female and married are factors)

model <- lm.ma(lwage ~ female + married + educ + exper + tenure, compute.deriv = TRUE, data = wage1)
summary(model)

## Compare coefficients from the linear model with the mean values from model averaging 
## for the non-factor predictors (but note in the summary that there is substantial variation)

apply(model$deriv,2,summary)

colMeans(model$deriv)

## Linear regression derivatives

coef(model.lm)[4:6]
}

\keyword{Regression}

