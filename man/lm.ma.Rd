\name{lm.ma}
\alias{lm.ma}
\alias{lm.ma.default}
\alias{lm.ma.formula}
\title{
Fitting Model Average Models
}
\description{
A function with an interface similar to \code{\link{lm}} that averages over a set of linear (in parameters) candidate models.
}
\usage{
lm.ma(...)

\method{lm.ma}{default}(y = NULL,
      X = NULL,
      X.eval = NULL,
      basis = c("auto","tensor","glp","additive"),
      compute.deriv = FALSE,
      deriv.order = 1,
      degree.min = 1,
      degree.max = NULL,
      lambda = 1e-02,
      segments.max = 3,
      knots = FALSE,
      S = 2,
      method = c("jma","mma"),
      ma.weights = NULL,
      basis.vec = NULL,
      bootstrap.ci = FALSE,
      B = 199,
      alpha = 0.05,
      weights = NULL,
      vc = TRUE,
      verbose = TRUE,
      tol = 1e-12,
      compute.anova = FALSE,
      ...)

\method{lm.ma}{formula}(formula,
      data = list(),
      y = NULL,
      X = NULL,
      X.eval = NULL,
      basis = c("auto","tensor","glp","additive"),
      compute.deriv = FALSE,
      deriv.order = 1,
      degree.min = 1,
      degree.max = NULL,
      lambda = 1e-02,
      segments.max = 3,
      knots = FALSE,
      S = 2,
      method = c("jma","mma"),
      ma.weights = NULL,
      basis.vec = NULL,
      bootstrap.ci = FALSE,
      B = 199,
      alpha = 0.05,
      weights = NULL,
      vc = TRUE,
      verbose = TRUE,
      tol = 1e-12,
      compute.anova = FALSE,
      ...)

}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit }

  \item{data}{ an optional data frame containing the variables in the model }

  \item{y}{a one dimensional vector of dependent data}

  \item{X}{a \eqn{p}-variate data frame of explanatory (training) data }

  \item{X.eval}{a \eqn{p}-variate data frame of points on which the regression will be estimated (evaluation data)}

  \item{basis}{ a character string indicating whether the generalized polynomial, additive or tensor product basis should be used}

  \item{compute.deriv}{a logical value indicating whether or not to compute derivatives}
  
  \item{deriv.order}{integer indicating order of derivative desired}

  \item{degree.min}{minimum value for the basis in each dimension}

  \item{degree.max}{maximum value for the basis in each dimension (defaults to \code{max(2,ceiling(S*log(NROW(X))/num.x))} where \code{num.x} is the number of numeric (i.e., non categorical) predictors and n the number of observations}

  \item{lambda}{smoothing parameter for weight function for the
  categorical predictors when \code{vc=TRUE}}
  
  \item{segments.max}{when \code{knots=TRUE} the maximum number of segments (i.e., interior knots plus 1) to use for the B-spline basis}
  
  \item{knots}{a logical value indicating whether to include interior knots or not (set equal to the degree)}

  \item{S}{the constant in the data-driven rule for determining degree.max}
  
  \item{method}{whether to use jackknife model averaging (jma, Hansen and Racine (2013)) or Mallows model averaging (mma, Hansen (2007))}

  \item{ma.weights}{a vector of model average weights obtained from a previous invocation (useful for bootstrapping etc.)}
  
  \item{basis.vec}{a vector (character) of bases for each candidate model selected by cross-validation}
  
  \item{bootstrap.ci}{a logical value indicating whether to bootstrap nonparametric confidence intervals or not}
  
  \item{B}{the number of bootstrap replications desired}
  
  \item{alpha}{value in (0,1) used to compute \eqn{1-\alpha\%} confidence intervals}
  
  \item{weights}{an optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector; if non-NULL, weighted least squares is used with weights \code{weights} (that is, minimizing \eqn{\sum_iw_ie_i^2}); otherwise ordinary least squares is used}

  \item{vc}{a logical value indicating whether to allow the categorical predictors to enter additively (only the intercept can shift) or to instead use a varying coefficient structure (all parameters can shift)}
  
  \item{verbose}{a logical value indicating whether to report detailed progress during computation}
  
  \item{tol}{a small constant used to ensure the forcing matrix for the quadratic program is not ill-conditioned}
  
  \item{compute.anova}{whether to conduct a quasi-snova bootstrap procedure to test for significance of each predictor}
  
  \item{...}{optional arguments to be passed}

}
\details{

Models for lm.ma are specified symbolically. A typical model has the form \code{response ~ terms} where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. Typical usages are
\preformatted{
    model <- lm.ma(y~x1+x2)
    model <- lm.ma(y~x1+x2,compute.deriv=TRUE)
    model <- lm.ma(y~x1+x2,bootstrap.ci=TRUE)  

    plot(model)
    plot(model,plot.data=TRUE)
    plot(model,plot.ci=TRUE,B=19)    
    plot(model,plot.data=TRUE,plot.ci=TRUE,B=199)
    plot(model,plot.deriv=TRUE)    
    plot(model,plot.deriv=TRUE,plot.ci=TRUE,B=399)        

    summary(model)
    fitted(model)
    coef(model)

    ## For generating predictions, create foo, a dataframe with named elements (important)
    ## for all predictors in the object model, then call predict, e.g.,
    foo <- data.frame(x1=c(1,2),x2=c(3,4))
    predict(model,newdata=foo)$fit

    ## If you want to see the degrees and knots selected by the procedure
    ## for the models that receive positive model average weights, try the
    ## following:
    model$DS[model$ma.weights>1e-05,]
}

Note that, unlike \code{lm} in which the formula interface specifies the functional form, in \code{lm.ma} the formula interface is strictly for listing the variables involved, while the bases will determine an appropriate model averaged functional form. Do not incorporate transformations, interactions and the like in the formula interface for \code{lm.ma} as these will most surely fail.

This function computes model averages over least squares candidate models whose predictors are generated by common basis functions (additive, multivariate polynomial, or tensor products). The candidate models increase in complexity from linear bases through higher order ones up to \code{degree.max}. The basis function (\code{basis="glp"}) uses a (multivariate) Taylor-type basis with Bernstein (as opposed to raw) polynomials. When \code{knots=TRUE} then interior knots are used and the Bernstein polynomials become B-spline bases and we are then averaging over regression spline models. For degrees higher than one, the multivariate polynomial basis includes interaction terms up to order degree minus one. Since we are averaging over models that are nonlinear in the predictors, derivatives will be vectors that potentially depend on the values of all predictors. An ad-hoc formula is used to determine the relationship between the largest (most complex) model, the sample size, and the number of predictors. This was chosen so that, as the sample size increases, we can approximate ever more complex functions while we must necessarily restrict the size of the largest model in small sample settings. Categorical predictors can enter additively and linearly (\code{vc=FALSE}) or in a parsimonious manner by exploiting recent developments in semiparametric varying coefficient models along the lines of Li, Ouyang, and Racine (2013). With the options \code{knots=TRUE} and \code{vc=TRUE}, we are averaging over varying-coefficient regression splines.

This approach frees the user from using either model assertion or selection methods and thereby attenuates bias arising from model misspecification. Simulations reveal that this approach is competitive with the best existing semi- and nonparametric approaches. Because it uses only least squares fits, it is also fairly computationally efficient.

}
\value{

\code{lm.ma} returns an object of class "lm.ma".

The function summary is used to obtain and print a summary of the results. The generic accessor functions \code{coef}, \code{fitted}, \code{predict}, \code{plot} and \code{residuals} extract various useful features of the value returned by lm.ma.

An object of class "lm.ma" is a list containing at least the following components:

  \item{fitted.values}{vector of fitted values}
  \item{deriv}{matrix of derivative vectors for each predictor}
  \item{ma.weights}{model average weights}
  \item{degree.max}{value of degree.max for each dimension (set by the default rule given above unless manually overridden)}
  \item{fitted.ci.l}{\eqn{\alpha/2} nonparametric confidence value vector for the vector of fitted/predicted values}
  \item{fitted.ci.u}{\eqn{1-\alpha/2} nonparametric confidence value vector for the vector fitted/predicted values} 
  \item{fitted.scale}{robust scale (mad) vector for the vector fitted/predicted values} 
  \item{deriv.ci.l}{\eqn{\alpha/2} nonparametric confidence value matrix for the matrix of derivatives}
  \item{deriv.ci.u}{\eqn{1-\alpha/2} nonparametric confidence value matrix for the matrix of derivatives}
  \item{deriv.scale}{robust scale matrix for the matrix of derivatives}
  \item{r.squared}{appropriate measure of goodness of fit (Doksum and Samarov (1995))}
  \item{residuals}{model residuals}

}
\references{
  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric estimation of global functionals and a measure of the explanatory power of covariates in regression,} The Annals of Statistics, 23 1443-1473.  
  
  Li, Q. and D. Ouyang and J.S. Racine (2013), \dQuote{Categorical Semiparametric Varying Coefficient Models,} Journal of Applied Econometrics, Volume 28, 551-579.
  
  Hansen, B. E. (2007), \dQuote{Least Squares Model Averaging,} Econometrica 75, 1175-1189.

  Hansen, B. E. & Racine, J. S. (2012), \dQuote{Jackknife Model Averaging,} Journal of Econometrics 167(1), 38-46.

  Racine, J.S. and D. Zhang and Q. Li (2017), \dQuote{Model Averaged Categorical Regression Splines.}
}
\author{
Jeffrey S. Racine
}
\note{
This code is in beta status until further notice - proceed accordingly.
}

\seealso{
\code{\link{lm}}, \code{\link[crs]{crs}}, \code{\link[np]{npreg}}
}
\examples{
#### Example 1 - simulated nonlinear one-predictor function

set.seed(42)
n <- 100
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + rnorm(n,sd=0.5*sd(dgp))

model.ma <- lm.ma(y~x)

summary(model.ma)

par(mfrow=c(1,2))
plot(model.ma,plot.data=TRUE,plot.ci=TRUE)
plot(model.ma,plot.data=TRUE,plot.ci=TRUE,plot.deriv=TRUE)
par(mfrow=c(1,1))

#### Example 2 - simulated multiplicative nonlinear two-predictor function

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

model.ma <- lm.ma(y~x1+x2)

summary(model.ma)

## Use the rgl package to render a 3D object (RGL is a 3D real-time rendering 
## system for R that supports OpenGL, among other formats).

z <- matrix(predict(model.ma,newdata=newdata)$fit,n.eval,n.eval)
num.colors <- 1000
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]
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

#### Example 3 - five predictor (two categorical) earnings function

data(wage1)
attach(wage1)

## Classical linear regression model (linear, additive, no interactions)

model.lm <- lm(lwage ~ female + married + educ + exper + tenure, 
               data = wage1)

## Murphy-Welch's favourite specification               
model.lm.mw <- lm(lwage ~ female + married + educ + exper + I(exper^2) 
                  + I(exper^3) + I(exper^4) + tenure)   

## Murphy-Welch's favourite specification with interactions in the intercepts              
model.lm.mwint <- lm(lwage ~ female + married + female:married + educ + exper 
                     + I(exper^2) + I(exper^3) + I(exper^4) + tenure)     

summary(model.lm)

## Compare with the model average estimator (female and married are factors)

model.ma <- lm.ma(lwage ~ female + married + educ + exper + tenure,
                  compute.deriv = TRUE)

summary(model.ma)

## Compare coefficients from the simple linear model with the (vector summary) values 
## from model averaging for the non-factor predictors

apply(coef(model.ma),2,summary)
coef(model.lm)[4:6]

## Compute parametric and model averaged marriage premiums for males and females
## at median values of remaining predictors

newdata.female.married <- data.frame(educ=round(median(educ)),
                                     exper=round(median(exper)),
                                     tenure=round(median(tenure)),
                                     female=factor("Female",levels=levels(female)),
                                     married=factor("Married",levels=levels(married)))

newdata.female.notmarried <- data.frame(educ=round(median(educ)),
                                        exper=round(median(exper)),
                                        tenure=round(median(tenure)),
                                        female=factor("Female",levels=levels(female)),
                                        married=factor("Notmarried",levels=levels(married)))
                                     
## Compute the so-called marriage premium - try three simple parametric 
## specifications (take your pick - is the premium +13%? +3%? -12%?) 

## Linear parametric
predict(model.lm,newdata=newdata.female.married)-
predict(model.lm,newdata=newdata.female.notmarried)

## Murphy-Welch parametric
predict(model.lm.mw,newdata=newdata.female.married)-
predict(model.lm.mw,newdata=newdata.female.notmarried)

## Murphy-Welch parametric augmented with a dummy interaction
predict(model.lm.mwint,newdata=newdata.female.married)-
predict(model.lm.mwint,newdata=newdata.female.notmarried)

## Model average
predict(model.ma,newdata=newdata.female.married)$fit-
predict(model.ma,newdata=newdata.female.notmarried)$fit

detach(wage1)

#### Example 4 - Canadian Current Population Survey earnings data

suppressPackageStartupMessages(require(np))
suppressPackageStartupMessages(require(crs))
data(cps71)
attach(cps71)
model.ma <- lm.ma(logwage~age)
plot(model.ma,plot.data=TRUE)
model.kernel <- npreg(logwage~age,regtype="ll",bwmethod="cv.aic")
lines(age,fitted(model.kernel),col=4,lty=4,lwd=2)
model.spline <- crs(logwage~age,cv.threshold=0)
lines(age,fitted(model.spline),col=3,lty=3,lwd=2)
legend("topleft",c("Model Average",
                   "Nonparametric Kernel",
                   "Nonparametric B-Spline"),
                   col=c(1,4,3),
                   lty=c(1,4,3),
                   lwd=c(1,2,2),
                   bty="n")
                   
summary(model.spline)
summary(model.kernel)
summary(model.ma)

detach(cps71)
}

\keyword{Regression}

