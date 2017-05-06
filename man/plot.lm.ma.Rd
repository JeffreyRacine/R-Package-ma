\name{plot.lm.ma}
\alias{plot.lm.ma}
\title{Plot an \code{lm.ma} Object}
\description{Plots a model average model and its derivatives.}
\usage{
\method{plot}{lm.ma}(x,
     plot.B = 99,
     plot.ci = FALSE,
     plot.data = FALSE,
     plot.deriv = FALSE,
     plot.num.eval = 250,
     plot.rug = FALSE,
     plot.xtrim = 0.005,
     ...)
}

\arguments{
  \item{x}{an object of type \code{lm.ma}}
  \item{plot.B}{number of bootstrap replications used to construct nonparametric confidence intervals}
  \item{plot.ci}{a logical value indicating whether to plot nonparametric confidence intervals or not}
  \item{plot.data}{a logical value indicating whether to plot the data or not}
  \item{plot.deriv}{a logical value indicating whether to compute derivatives or not}
  \item{plot.num.eval}{number of evaluation points}
  \item{plot.rug}{a logical value indicating whether to plot the data with a rug or not}
  \item{plot.xtrim}{trimming parameter used to exclude tail values for the predictors that can obscure main features in the plot (trims the proportion \code{plot.xtrim} from each tail)}
  \item{\dots}{optional arguments to be passed to \code{plot}}
}
\details{
This function plots an object returned by \code{lm.ma}. Typical usages are
\preformatted{
    plot(model)
    plot(model,plot.data=TRUE)
    plot(model,plot.ci=TRUE,plot.B=99)    
    plot(model,plot.data=TRUE,plot.ci=TRUE,plot.B=199)
    plot(model,plot.deriv=TRUE)    
    plot(model,plot.deriv=TRUE,plot.ci=TRUE,plot.B=399)        
  }
}
\value{
None.
}
\references{
Racine, J.S. and D. Zhang and Q. Li (2017), \dQuote{Model Averaged Categorical Regression Splines.}
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
