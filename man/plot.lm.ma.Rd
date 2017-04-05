\name{plot.lm.ma}
\alias{plot.lm.ma}
\title{Plot an \code{lm.ma} object}
\description{Plots a model average model and its derivatives.}
\usage{
\method{plot}{lm.ma}(x,
     B = 99,
     plot.deriv = FALSE,
     plot.ci = FALSE,
     plot.data = FALSE,
     plot.num.eval = 250,
     plot.xtrim = 0.001,
     ...)
}

\arguments{
  \item{x}{an object of type \code{lm.ma}}
  \item{B}{number of bootstrap replications used to construct nonparametric confidence intervals}
  \item{plot.deriv}{a logical value indicating whether to compute derivatives or not}
  \item{plot.ci}{a logical value indicating whether to plot nonparametric confidence intervals or not}
  \item{plot.data}{a logical value indicating whether to plot the data or not}
  \item{plot.num.eval}{number of evaluation points}
  \item{plot.xtrim}{small trimming parameter to exclude a few outliers that can obscure main features in the plot}
  \item{\dots}{optional arguments to be passed to \code{plot}}
}
\details{
This function plots an object returned by \code{lm.ma}.
}
\value{
None.
}
\references{
Racine, J.S. and D. Zhang and Q. Li (2017), \dQuote{Model Averaged Categorical Regression Splines.}
}
None.
}
\author{
Jeffrey S. Racine
}
\examples{
data(cps71)
model <- lm.ma(logwage~age,data=cps71)
plot(model,plot.data=TRUE,plot.ci=TRUE)
}

\keyword{Regression}
