\name{lm.ma}
\alias{lm.ma}
\alias{lm.ma.default}
\alias{lm.ma.formula}
\title{Fitting Model Average Models}
\description{
A function with an interface similar to \code{\link{lm}} that averages over a set of linear (in parameters) candidate models.
}
\usage{
lm.ma(...)

\method{lm.ma}{default}(y = NULL,
      X = NULL,
      X.eval = NULL,
      all.combinations = TRUE,
      alpha = 0.05,
      auto.basis = c("tensor","taylor","additive"),
      auto.reduce = TRUE,
      B = 99,
      basis.vec = NULL,
      basis = c("auto","tensor","taylor","additive"),
      boot.ci = FALSE,
      compute.anova = FALSE,
      compute.anova.boot = FALSE,
      compute.anova.index = NULL,
      compute.deriv = FALSE,
      compute.mean = TRUE,
      degree.by = 2,
      degree.max = NULL,
      degree.min = 0,
      deriv.index = NULL,
      deriv.order = 1,
      DKL.mat = NULL,
      eps.lambda = 1e-04,
      knots = FALSE,
      lambda.S = 2,
      lambda.max = 1,
      lambda.num.max = NULL,
      ma.weights = NULL,
      ma.weights.cutoff = 1e-04,
      max.dim.candidate.models = 5000,
      max.num.candidate.models = 2500,
      method = c("jma","mma"),
      parallel = FALSE,
      parallel.cores = NULL,
      rank.vec = NULL,
      restrict.sum.ma.weights = TRUE,
      rng.seed = 42,
      S = 1,
      segments.by = 2,
      segments.max = 3,
      segments.min = 1,
      singular.ok = TRUE,
      trace = FALSE,
      vc = TRUE,
      verbose = TRUE,
      weights = NULL,
      ...)

\method{lm.ma}{formula}(formula,
      data = list(),
      y = NULL,
      X = NULL,
      X.eval = NULL,
      all.combinations = TRUE,
      alpha = 0.05,
      auto.basis = c("tensor","taylor","additive"),
      auto.reduce = TRUE,
      B = 99,
      basis.vec = NULL,
      basis = c("auto","tensor","taylor","additive"),
      boot.ci = FALSE,
      compute.anova = FALSE,
      compute.anova.boot = FALSE,
      compute.anova.index = NULL,
      compute.deriv = FALSE,
      compute.mean = TRUE,
      degree.by = 2,
      degree.max = NULL,
      degree.min = 0,
      deriv.index = NULL,
      deriv.order = 1,
      DKL.mat = NULL,
      eps.lambda = 1e-04,
      knots = FALSE,
      lambda.S = 2,
      lambda.max = 1,
      lambda.num.max = NULL,
      ma.weights = NULL,
      ma.weights.cutoff = 1e-04,
      max.dim.candidate.models = 5000,
      max.num.candidate.models = 2500,
      method = c("jma","mma"),
      parallel = FALSE,
      parallel.cores = NULL,
      rank.vec = NULL,
      restrict.sum.ma.weights = TRUE,
      rng.seed = 42,
      S = 1,
      segments.by = 2,
      segments.max = 3,
      segments.min = 1,
      singular.ok = TRUE,
      trace = FALSE,
      vc = TRUE,
      verbose = TRUE,
      weights = NULL,
      ...)

}

\arguments{
  \item{formula}{a symbolic description of the model to be fit }

  \item{data}{an optional data frame containing the variables in the model }

  \item{y}{a one dimensional vector of dependent data}

  \item{X}{a \eqn{p}-variate data frame of explanatory (training) data }

  \item{X.eval}{a \eqn{p}-variate data frame of points on which the regression will be estimated (evaluation data)}

  \item{all.combinations}{a logical value indicating whether or not to attempt all combinations of degrees, segments, knots, and lambda values (if \code{all.combinations=FALSE} only candidate models with the same degree in all dimensions are considered, no interior knots are considered, while the minimum number of segments and smallest value for lambda are used)}

  \item{alpha}{a value in (0,1) used to compute \eqn{1-\alpha\%} confidence intervals}

  \item{auto.basis}{which (subset possible) bases to use when \code{basis="auto"}}
  
  \item{auto.reduce}{a logical value indicating whether or not to use some crude heuristics to reduce the number of candidate models if the number of candidate models exceeds \code{max.num.candidate.models}}
  
  \item{B}{the number of bootstrap replications desired}

  \item{basis.vec}{a vector (character) of bases for each candidate model}

  \item{basis}{a character string indicating whether the generalized Taylor polynomial, additive or tensor product basis should be used (if \code{basis="auto"} then for each candidate model the most appropriate basis is determined via cross-validation)}
  
  \item{boot.ci}{a logical value indicating whether or not to construct nonparametric bootstrap confidence intervals}

  \item{compute.anova}{a logical value indicating whether or not to conduct an anova-based procedure to test for predictor significance}

  \item{compute.anova.boot}{a logical value indicating whether or not the test for predictor significance uses asymptotic or bootstrapped P-values}

  \item{compute.anova.index}{an optional vector of indices indicating which predictor(s) are to be tested (default is all predictors)}

  \item{compute.deriv}{a logical value indicating whether or not to compute derivatives}

  \item{compute.mean}{a logical value indicating whether or not to
    compute the conditional mean}

  \item{degree.by}{increment in degree sequence (if \code{degree.min=0} sequence will include degree 0 then start at 1 in increments of \code{degree.by})}

  \item{degree.max}{the maximum value for the basis degree in each dimension (the value defaults to \code{max(2,ceiling(log(n)-S*log(1+k)))} where \code{k} is the number of numeric predictors and \code{n} the number of observations}

  \item{degree.min}{the minimum value for the basis degree in each dimension}

  \item{deriv.index}{an optional vector of indices indicating which predictor(s) derivative is computed}
  
  \item{deriv.order}{an integer indicating the order of derivative desired (1,2,...)}

  \item{DKL.mat}{a matrix with degree, knots, and lambda values (if \code{vc=TRUE}) that could optionally be passed to the basis routines}

  \item{eps.lambda}{a small positive constant for the start of the sequence of smoothing parameters used in the weight function for the categorical predictors when \code{vc=TRUE}}

  \item{knots}{a logical value indicating whether or not to include interior knots}

  \item{lambda.S}{the constant in the data-driven rule for determining \code{lambda.num.max}}

  \item{lambda.max}{largest value (<= 1) of the smoothing parameters used in the weight function for the categorical predictors when \code{vc=TRUE}}

  \item{lambda.num.max}{the maximum value for the smoothing parameter grid in each dimension (defaults to \code{max(2,ceiling(log(n)-lambda.S*log(1+p)))} where \code{p} is the number of categorical predictors and \code{n} the number of observations}
  
  \item{ma.weights}{a vector of model average weights obtained from a previous invocation (useful for bootstrapping etc.)}

  \item{ma.weights.cutoff}{a small number below which a model weight is deemed to be essentially zero}

  \item{max.dim.candidate.models}{an arbitrary upper bound on the maximum dimension of candidate models permitted}
  
  \item{max.num.candidate.models}{an arbitrary upper bound on the maximum number of candidate models permitted}
  
  \item{method}{a character string indicating whether to use jackknife model averaging (\code{"jma"}, Hansen and Racine (2013)) or Mallows model averaging (\code{"mma"}, Hansen (2007) - both are frequentist model average criterion)}

  \item{parallel}{a logical value indicating whether or not to run certain routines in parallel}
  
  \item{parallel.cores}{a positive integer indicating the number of cores desired when \code{parallel=TRUE} (when \code{parallel=FALSE} defaults to the number of available cores)}
  
  \item{rank.vec}{a vector of ranks for each candidate model}

  \item{restrict.sum.ma.weights}{a logical value indicating whether or not to restrict the sum of the model average weights to one when solving the quadratic program (they are normalized afterwards when \code{restrict.sum.ma.weights=FALSE})}
  
  \item{rng.seed}{an integer used to seed R's random number generator - this is to ensure replicability when bootstrapping}

  \item{S}{the constant in the data-driven rule for determining \code{degree.max}}

  \item{segments.by}{increment in segments sequence when \code{knots=TRUE}}

  \item{segments.min}{the minimum number of segments when \code{knots=TRUE} (i.e., number of knots minus 1 - there always exist two knots, the endpoints) to allow for the B-spline basis}

  \item{segments.max}{the maximum number of segments when \code{knots=TRUE} (i.e., number of knots minus 1 - there always exist two knots, the endpoints) to allow for the B-spline basis}

  \item{singular.ok}{if \sQuote{FALSE} (the default in S but not in R) a singular fit is an error}
  
  \item{trace}{a logical value indicating whether or not to issue a detailed progress report via \code{warning}}
  
  \item{vc}{a logical value indicating whether to allow the categorical predictors to enter additively (only the intercept can shift) or to instead use a varying coefficient structure (all parameters can shift)}

  \item{verbose}{a logical value indicating whether to report detailed progress during computation (\code{warnings()} are not issued when \code{verbose=FALSE})}

  \item{weights}{an optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector; if non-NULL, weighted least squares is used with weights \code{weights} (that is, minimizing \eqn{\sum_iw_ie_i^2}); otherwise ordinary least squares is used}

  \item{...}{optional arguments to be passed}

}
\details{

Models for lm.ma are specified symbolically. A typical model has the form \code{response ~ terms} where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. Typical usages are
\preformatted{
    model <- lm.ma(y~x1+x2)
    model <- lm.ma(y~x1+x2,compute.deriv=TRUE)
    model <- lm.ma(y~x1+x2,boot.ci=TRUE)  
    model <- lm.ma(y~x1+x2,compute.anova=TRUE,compute.anova.boot=TRUE,degree.min=1)
    model <- lm.ma(y~x1+x2,parallel=TRUE)
    model <- lm.ma(y~x1+x2,parallel=TRUE,parallel.cores=2)    
    model <- lm.ma(y~x+z,lambda.S=3)
    model <- lm.ma(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                   vc=FALSE,
                   degree.by=1,
                   degree.max=5,
                   basis="additive",
                   all.combinations=FALSE)

    plot(model)
    plot(model,plot.data=TRUE)
    plot(model,plot.ci=TRUE,plot.B=199)    
    plot(model,plot.data=TRUE,plot.ci=TRUE,plot.B=199)
    plot(model,plot.deriv=TRUE)    
    plot(model,plot.deriv=TRUE,plot.ci=TRUE,plot.B=399)        

    summary(model)
    fitted(model)
    coef(model)

    ## For generating predictions, create foo, a data frame with named
    ## elements (important) for all predictors in the object model,
    ## then call predict, e.g.,

    foo <- data.frame(x1=c(1,2),x2=c(3,4))
    predict(model,newdata=foo)

    ## If you want to see the degrees, number of segments, and smoothing
    ## parameters for the categorical predictors (vc=TRUE) selected by
    ## the procedure for the models that receive positive model average
    ## weights, try the following:

    model$DKL.mat[model$ma.weights>model$ma.weights.cutoff,]
}

Note that, unlike \code{lm} in which the formula interface specifies functional form, in \code{lm.ma} the formula interface is strictly for listing the variables involved and the procedure will determine an appropriate model averaged functional form. Do not incorporate transformations, interactions and the like in the formula interface for \code{lm.ma} as these will most surely fail.

This function computes a model that is the weighted average of a set of least squares candidate models whose predictors are generated by common basis functions (additive, generalized Taylor polynomial, or tensor products). The candidate models increase in complexity from linear bases (if \code{degree.min=1}) through higher order ones up to \code{degree.max}. All bases are of the Bernstein polynomial class, as opposed to raw polynomials, and allow for differing degrees across multivariate predictors. When \code{knots=TRUE}, interior knots are used and the Bernstein polynomials become B-spline bases and we are then averaging over regression spline models. When the number of numeric predictors is two or more, the generalized Taylor polynomial includes interaction terms up to order degree minus one. Since we are averaging over models that are nonlinear in the predictors, derivatives will be vectors that potentially depend on the values of every predictor. An ad-hoc formula is used to determine the relationship between the largest (most complex) model, the sample size, and the number of predictors. This ad-hoc rule was set so that, as the sample size increases, we can approximate ever more complex functions while necessarily restricting the size of the largest model in small sample settings. Categorical predictors can enter additively and linearly (if \code{vc=FALSE}) or in a parsimonious manner by exploiting recent developments in semiparametric varying coefficient models along the lines of Li, Ouyang, and Racine (2013). With the options \code{knots=TRUE} and \code{vc=TRUE}, we are averaging over varying-coefficient regression splines.

This approach frees the user from using either model assertion or selection methods and thereby attenuates bias arising from model misspecification. Simulations reveal that this approach is competitive with some semi- and nonparametric approaches. Because it uses only least squares fits, it can be more computationally efficient than its nonparametric counterparts.

}
\value{

\code{lm.ma} returns an object of class "lm.ma".

The function summary is used to obtain and print a summary of the results. The generic accessor functions \code{coef}, \code{fitted}, \code{predict}, \code{plot} (see \code{?plot.lm} for details) and \code{residuals} extract various useful features of the value returned by lm.ma.

An object of class "lm.ma" is a list containing at least the following components:
  \item{degree.max}{value of degree.max for each dimension (set by an ad-hoc rule unless manually overridden)}
  \item{deriv.ci.l}{\eqn{\alpha/2} nonparametric confidence value matrix for the matrix of derivatives}
  \item{deriv.ci.u}{\eqn{1-\alpha/2} nonparametric confidence value matrix for the matrix of derivatives}
  \item{deriv.scale}{robust scale (mad) matrix for the matrix of derivatives}
  \item{deriv}{matrix of derivative vectors for each predictor}
  \item{fitted.ci.l}{\eqn{\alpha/2} nonparametric confidence value vector for the vector of fitted/predicted values}
  \item{fitted.ci.u}{\eqn{1-\alpha/2} nonparametric confidence value vector for the vector fitted/predicted values}
  \item{fitted.scale}{robust scale (mad) vector for the vector of fitted/predicted values}
  \item{fitted.values}{vector of fitted values}
  \item{ma.weights}{model average weights}
  \item{r.squared}{appropriate measure of goodness of fit (Doksum and Samarov (1995))}
  \item{residuals}{model residuals}

}
\references{
  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric Estimation of Global Functionals and a Measure of the Explanatory Power of Covariates in Regression,} The Annals of Statistics, 23 1443-1473.  
  
  Li, Q. and D. Ouyang and J.S. Racine (2013), \dQuote{Categorical Semiparametric Varying Coefficient Models,} Journal of Applied Econometrics, Volume 28, 551-579.
  
  Hansen, B. E. (2007), \dQuote{Least Squares Model Averaging,} Econometrica 75, 1175-1189.

  Hansen, B. E. & Racine, J. S. (2012), \dQuote{Jackknife Model Averaging,} Journal of Econometrics 167(1), 38-46.

  Racine, J.S. and Q. Li and L. Zheng (2018), \dQuote{Optimal Model Averaging of Mixed-Data Kernel-Weighted Spline Regressions.}
}
\author{
Jeffrey S. Racine
}
\note{
This code is in beta status until further notice - proceed accordingly.

Note that the purpose of this package is to attenuate bias arising from model misspecification in situations where model uncertainty is present and you are concerned about its impact on any subsequent inference and prediction. This package is best suited to situations involving a manageable number of predictors (i.e., a handful or two at most) and a sufficient number of observations so that nonlinearities can reasonably be uncovered. If your objective is to include all possible measured predictors (i.e., the kitchen sink approach) and conduct variable selection (i.e., attempt to determine which variables enter linearly), this package is not for you; see instead the R packages \pkg{BMA}, \pkg{lars}, or the function \link[MASS]{stepAIC} in the \pkg{MASS} package (with \code{degree.max=1} the defaults would only allow for at most eleven numeric predictors, i.e., \eqn{2^{11}} combinations of degrees 0 and 1). To get around this limitation that arises by attempting to consider a range of degree, segment, and smoothing parameter values for each dimension (the number of combinations can quickly get far too large), the option \code{all.combinations=FALSE} can be invoked. This restricts the number of candidate models by holding the degree, segment, and smoothing parameters to be the same for each dimension which can reduce the number of models to just a handful at most. Using \code{basis="additive"} further restricts the rank of each candidate model, while \code{vc=FALSE} can reduce execution time in the presence of categorical predictors (see the example in \bold{Details} above).

The number of candidate models may grow unreasonably large (say 2,500 or more) if multiple predictors are present. Some heuristics are therefore necessary in order to corral the number of candidate models (and the maximum basis dimension). However, no default setting can be ideal for all data generating processes and you may wish to intervene. If you wish to reduce the number of candidate models used, there are a number of ways of accomplishing this. In particular, you might want to i) increase \code{S}, ii) increase \code{lambda.S} (if categorical predictors are present and \code{vc=TRUE}), iii) set and restrict \code{degree.max}, iv) set and restrict \code{lambda.num.max} if categorical predictors are present, v) reduce \code{segments.max} (if \code{knots=TRUE}),  vi) set \code{all.combinations=FALSE}, vi) directly modify \code{max.dim.candidate.models} and/or \code{max.num.candidate.models}, or perhaps instead consider a semiparametric model (\code{basis="additive"} and \code{vc=FALSE} produces semiparametric additive candidate models - see the example in \code{?india} for an illustration). When building the final model each candidate model must be constructed and evaluated. However, after solving for the model average weights, a number of candidate models may be assigned essentially zero weight. Subsequently, only the non-zero weight models need be evaluated (e.g. when constructing derivative estimates, predictions, confidence intervals and the like).

When \code{compute.anova.boot=TRUE}, the option \code{compute.anova} uses a bootstrap procedure that requires re-computation of the model average model for each bootstrap replication. With one or two predictors and \code{compute.anova.boot=TRUE} the procedure may be fairly fast, but as the model complexity increases the procedure will require some patience.

The option \code{compute.anova=TRUE} cannot be used in the presence of one or more factors and exactly one numeric predictor since there is no numeric predictor present when testing for significance for the one numeric predictor.

Note that the option \code{compute.anova=TRUE} (not default) will warn immediately when \code{degree.min=0} (default) and rest to \code{degree.min=1}. The reason for this is because irrelevant predictors can be automatically removed without the need for pre-testing if the procedure selects the degree for any predictor to be 0 - in such cases the restricted and unrestricted models may coincide and the test is degenerate. The same holds for smoothing parameter values with \code{vc=TRUE} in the presence of categorical predictors (when \code{lambda=1} irrelevant categorical predictors are automatically removed without the need for pretesting, so we need to rule this case out when conducting hypothesis tests).

Averaging over models with ill-conditioned bases is not advised. Pay attention to the warning \dQuote{Dimension basis is ill-conditioned - reduce \code{degree.max}} should it arise and reduce \code{degree.max} until this no longer is the case.

If you wish to plot the object with the option \code{plot.ci=TRUE}, it is not necessary to use the option \code{boot.ci=TRUE} in the call to \code{lm.ma()} (this will simply add to overhead)

Note that \code{predict.ma} produces a vector of predictions or a list of predictions, confidence bounds, derivative matrices and their confidence bounds.

See the examples contained in \code{demo(package="ma")} for illustrative demonstrations with real and simulated data (e.g., \code{demo(cps71,package="ma")}).

}

\seealso{
\code{\link{lm}}, \code{\link[crs]{crs}}, \code{\link[np]{npreg}}
}
\examples{
options(warn=-1)
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

## Compare with a semiparametric additive model average estimator
## (female and married are factors)

model.ma <- lm.ma(lwage ~ female + married + educ + exper + tenure,
                  compute.deriv = TRUE,
                  basis = "additive",
                  degree.by = 1,
                  vc = FALSE)

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

