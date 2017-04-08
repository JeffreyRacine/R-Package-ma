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
      alpha = 0.05,
      B = 99,
      B.scale = 0,
      basis.vec = NULL,
      basis = c("auto","tensor","taylor","additive"),
      bootstrap.ci = FALSE,
      compute.sigtest.anova = FALSE,
      compute.sigtest.deriv = FALSE,
      compute.deriv = FALSE,
      degree.max = NULL,
      degree.min = 1,
      deriv.index = NULL,
      deriv.order = 1,
      K.mat = NULL,
      knots = FALSE,
      lambda = 1e-02,
      ma.weights = NULL,
      method = c("jma","mma"),
      rank.vec = NULL,
      restrict.sum.ma.weights = FALSE,
      S = 2,
      segments.max = 3,
      tol = 1e-12,
      vc = TRUE,
      verbose = TRUE,
      weights = NULL,
      ...)

\method{lm.ma}{formula}(formula,
      data = list(),
      y = NULL,
      X = NULL,
      X.eval = NULL,
      alpha = 0.05,
      B = 99,
      B.scale = 0,
      basis.vec = NULL,
      basis = c("auto","tensor","taylor","additive"),
      bootstrap.ci = FALSE,
      compute.sigtest.anova = FALSE,
      compute.sigtest.deriv = FALSE,
      compute.deriv = FALSE,
      degree.max = NULL,
      degree.min = 1,
      deriv.index = NULL,
      deriv.order = 1,
      K.mat = NULL,
      knots = FALSE,
      lambda = 1e-02,
      ma.weights = NULL,
      method = c("jma","mma"),
      rank.vec = NULL,
      restrict.sum.ma.weights = FALSE,
      S = 2,
      segments.max = 3,
      tol = 1e-12,
      vc = TRUE,
      verbose = TRUE,
      weights = NULL,
      ...)

}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit }

  \item{data}{ an optional data frame containing the variables in the model }

  \item{y}{a one dimensional vector of dependent data}

  \item{X}{a \eqn{p}-variate data frame of explanatory (training) data }

  \item{X.eval}{a \eqn{p}-variate data frame of points on which the regression will be estimated (evaluation data)}

  \item{alpha}{value in (0,1) used to compute \eqn{1-\alpha\%} confidence intervals}

  \item{B}{the number of bootstrap replications desired}

  \item{B.scale}{the number of nested bootstrap replications desired to obtain a robust estimate of derivative scale when conducting hypothesis tests (if set to zero nested bootstrapping is not performed)}

  \item{basis.vec}{a vector (character) of bases for each candidate model selected by cross-validation}

  \item{basis}{a character string indicating whether the generalized polynomial, additive or tensor product basis should be used}

  \item{bootstrap.ci}{a logical value indicating whether to bootstrap nonparametric confidence intervals or not}

  \item{compute.sigtest.anova}{whether to conduct an anova-based bootstrap procedure to test for significance of each predictor}

  \item{compute.sigtest.deriv}{whether to conduct a derivative-based bootstrap procedure to test for significance of each predictor}

  \item{compute.deriv}{a logical value indicating whether or not to compute derivatives}

  \item{degree.max}{maximum value for the basis in each dimension (defaults to \code{max(2,S*log(n)/k} where \code{k} is the number of numeric (i.e., non categorical) predictors and \code{n} the number of observations}

  \item{degree.min}{minimum value for the basis in each dimension}


  \item{deriv.index}{optional vector of indices indicating as to which predictor(s) the derivative computation is to be performed}
  
  \item{deriv.order}{integer indicating order of derivative desired}

  \item{K.mat}{matrix with degree and knots that could be passed to the basis routines}

  \item{knots}{a logical value indicating whether to include interior knots or not (set equal to the degree)}

  \item{lambda}{smoothing parameter for weight function for the categorical predictors when \code{vc=TRUE}}

  \item{ma.weights}{a vector of model average weights obtained from a previous invocation (useful for bootstrapping etc.)}

  \item{method}{whether to use jackknife model averaging (jma, Hansen and Racine (2013)) or Mallows model averaging (mma, Hansen (2007))}

  \item{rank.vec}{a vector of ranks for each candidate model}

  \item{restrict.sum.ma.weights}{a logical value indicating whether or not to restrict the sum of the model average weights to one when solving the quadratic program (they are normalized afterwards when \code{restrict.sum.ma.weights=FALSE})}
  
  \item{S}{the constant in the data-driven rule for determining degree.max}

  \item{segments.max}{when \code{knots=TRUE} the maximum number of segments (i.e., interior knots plus 1) to use for the B-spline basis}

  \item{tol}{a small constant used to ensure the forcing matrix for the quadratic program is not ill-conditioned}

  \item{vc}{a logical value indicating whether to allow the categorical predictors to enter additively (only the intercept can shift) or to instead use a varying coefficient structure (all parameters can shift)}

  \item{verbose}{a logical value indicating whether to report detailed progress during computation}

  \item{weights}{an optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector; if non-NULL, weighted least squares is used with weights \code{weights} (that is, minimizing \eqn{\sum_iw_ie_i^2}); otherwise ordinary least squares is used}

  \item{...}{optional arguments to be passed}

}
\details{

Models for lm.ma are specified symbolically. A typical model has the form \code{response ~ terms} where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. Typical usages are
\preformatted{
    model <- lm.ma(y~x1+x2)
    model <- lm.ma(y~x1+x2,compute.deriv=TRUE)
    model <- lm.ma(y~x1+x2,bootstrap.ci=TRUE)  
    model <- lm.ma(y~x1+x2,compute.sigtest.anova=TRUE)
    model <- lm.ma(y~x1+x2,compute.sigtest.deriv=TRUE)      

    plot(model)
    plot(model,plot.data=TRUE)
    plot(model,plot.ci=TRUE,B=19)    
    plot(model,plot.data=TRUE,plot.ci=TRUE,B=199)
    plot(model,plot.deriv=TRUE)    
    plot(model,plot.deriv=TRUE,plot.ci=TRUE,B=399)        

    summary(model)
    fitted(model)
    coef(model)

    ## For generating predictions, create foo, a dataframe with named
    ## elements (important) for all predictors in the object model,
    ## then call predict, e.g.,

    foo <- data.frame(x1=c(1,2),x2=c(3,4))
    predict(model,newdata=foo)

    ## If you want to see the degrees and knots selected by the procedure
    ## for the models that receive positive model average weights, try the
    ## following:

    model$DS[model$ma.weights>1e-05,]
}

Note that, unlike \code{lm} in which the formula interface specifies the functional form, in \code{lm.ma} the formula interface is strictly for listing the variables involved, while the bases will determine an appropriate model averaged functional form. Do not incorporate transformations, interactions and the like in the formula interface for \code{lm.ma} as these will most surely fail.

This function computes model averages over least squares candidate models whose predictors are generated by common basis functions (additive, multivariate polynomial, or tensor products). The candidate models increase in complexity from linear bases through higher order ones up to \code{degree.max}. All bases are of the Bernstein polynomial class, as opposed to raw polynomials, and allow for differing degrees across multivariate predictors. When \code{knots=TRUE}, interior knots are used and the Bernstein polynomials become B-spline bases and we are then averaging over regression spline models. When the number of numeric predictors is two or more, the multivariate polynomial includes interaction terms up to order degree minus one. Since we are averaging over models that are nonlinear in the predictors, derivatives will be vectors that potentially depend on the values of every predictor. An ad-hoc formula is used to determine the relationship between the largest (most complex) model, the sample size, and the number of predictors. This ad-hoc rule was set so that, as the sample size increases, we can approximate ever more complex functions while necessarily restricting the size of the largest model in small sample settings. Categorical predictors can enter additively and linearly (\code{vc=FALSE}) or in a parsimonious manner by exploiting recent developments in semiparametric varying coefficient models along the lines of Li, Ouyang, and Racine (2013). With the options \code{knots=TRUE} and \code{vc=TRUE}, we are averaging over varying-coefficient regression splines.

This approach frees the user from using either model assertion or selection methods and thereby attenuates bias arising from model misspecification. Simulations reveal that this approach is competitive with some semi- and nonparametric approaches. Because it uses only least squares fits, it can be more computationally efficient than its nonparametric counterparts.

}
\value{

\code{lm.ma} returns an object of class "lm.ma".

The function summary is used to obtain and print a summary of the results. The generic accessor functions \code{coef}, \code{fitted}, \code{predict}, \code{plot} (see \code{plot.lm} for details) and \code{residuals} extract various useful features of the value returned by lm.ma.

An object of class "lm.ma" is a list containing at least the following components:
  \item{degree.max}{value of degree.max for each dimension (set by an ad-hoc rule unless manually overridden)}
  \item{deriv.ci.l}{\eqn{\alpha/2} nonparametric confidence value matrix for the matrix of derivatives}
  \item{deriv.ci.u}{\eqn{1-\alpha/2} nonparametric confidence value matrix for the matrix of derivatives}
  \item{deriv.scale}{robust scale matrix for the matrix of derivatives}
  \item{deriv}{matrix of derivative vectors for each predictor}
  \item{fitted.ci.l}{\eqn{\alpha/2} nonparametric confidence value vector for the vector of fitted/predicted values}
  \item{fitted.ci.u}{\eqn{1-\alpha/2} nonparametric confidence value vector for the vector fitted/predicted values}
  \item{fitted.scale}{robust scale (mad) vector for the vector fitted/predicted values}
  \item{fitted.values}{vector of fitted values}
  \item{ma.weights}{model average weights}
  \item{r.squared}{appropriate measure of goodness of fit (Doksum and Samarov (1995))}
  \item{residuals}{model residuals}

Note that \code{predict.ma} produces a vector of predictions or a list of predictions, confidence bounds, derivative matrices and their confidence bounds.
  
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

The options \code{compute.sigtest.anova} and \code{compute.sigtest.deriv} use a bootstrap procedure that requires re-computation of the model average model for each bootstrap replication. With one or two predictors it will be fairly fast, but as the model complexity increases it may take some patience. Virtually all of the approaches required herein are parallel in nature and could be run on multiple cores. If there were enough demand for this feature I might be encouraged to add a parallel option.

The option \code{compute.sigtest.anova} cannot be used in the presence of one or more factors and exactly one numeric predictor since there is no numeric predictor present when testing for significance for the one numeric predictor.

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

## Note that the following calls to plot() use the option
## plot.ci=TRUE which then invokes a bootstrap procedure. The
## plots may take a few seconds to appear due to this additional
## computation (if you remove this option the plots will appear
## sooner).

par(mfrow=c(1,2))
plot(model.ma,plot.data=TRUE,plot.ci=TRUE)
plot(model.ma,plot.data=TRUE,plot.ci=TRUE,plot.deriv=TRUE)
par(mfrow=c(1,1))

#### Example 2 - five predictor (two categorical) earnings function

data(wage1)
attach(wage1)

## Classical linear regression model (linear, additive, no interactions)

model.lm <- lm(lwage ~ female + married + educ + exper + tenure)

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

#### Example 3 - Canadian Current Population Survey earnings data

## We compute two nonparametric estimators to compare with the 
## model averaging approach. 

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

#### Example 5 - simulated multiplicative nonlinear two-predictor function

suppressPackageStartupMessages(require(rgl))

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

z <- matrix(predict(model.ma,newdata=newdata),n.eval,n.eval)
num.colors <- 1000
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]
par(ask=TRUE)
readline(prompt = "Hit <Return> to see next plot:")
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

## Note - if you click on the rgl window you can rotate the estimate 
## by dragging the object, zoom in and out etc.
}

\keyword{Regression}

