\name{india}
\alias{india}
\docType{data}
\title{
Childhood Malnutrition in India
}
\description{
Demographic and Health Survey data on childhood nutrition in India.
}
\usage{data(india)}
\format{
  A data frame with 37623 observations on the following 21 variables.
  \describe{
    \item{\code{cheight}}{child's height (centimeters); a numeric vector}
    \item{\code{cage}}{child's age (months); a numeric vector}
    \item{\code{breastfeeding}}{duration of breastfeeding (months); a numeric vector}
    \item{\code{csex}}{child's sex; a factor with levels \code{male} \code{female}}
    \item{\code{ctwin}}{whether or not child is a twin; a factor with levels \code{single birth} \code{twin}}
    \item{\code{cbirthorder}}{birth order of the child; a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
    \item{\code{mbmi}}{mother's BMI (kilograms per meter squared); a numeric vector}
    \item{\code{mage}}{mother's age (years); a numeric vector}
    \item{\code{medu}}{mother's years of education; a numeric vector}
    \item{\code{edupartner}}{father's years of education; a numeric vector}
    \item{\code{munemployed}}{mother's employment status; a factor variable with levels \code{unemployed} \code{employed}}
    \item{\code{mreligion}}{mother's religion; a factor variable with levels \code{christian} \code{hindu} \code{muslim} \code{other} \code{sikh}}
    \item{\code{mresidence}}{mother's residential classification; a factor with levels \code{urban} \code{rural}}
    \item{\code{wealth}}{mother's relative wealth; a factor with levels \code{poorest} \code{poorer} \code{middle} \code{richer} \code{richest}}
    \item{\code{electricity}}{electricity access; a factor with levels \code{no} \code{yes}}
    \item{\code{radio}}{radio ownership; a factor with levels \code{no} \code{yes}}
    \item{\code{television}}{television ownership; a factor with levels \code{no} \code{yes}}
    \item{\code{refrigerator}}{refrigerator ownership; a factor with levels \code{no} \code{yes}}
    \item{\code{bicycle}}{bicycle ownership; a factor with levels \code{no} \code{yes}}
    \item{\code{motorcycle}}{motorcycle ownership; a factor with levels \code{no} \code{yes}}
    \item{\code{car}}{car ownership; a factor with levels \code{no} \code{yes}}
  }
}
\source{
http://www.econ.uiuc.edu/~roger/research/bandaids/india.Rda
}
\references{
Koenker, R. (2011), "Additive models for quantile regression: Model selection and confidence bandaids," Brazilian Journal of Probability and Statistics 25(3), pp. 239-262.
}
\examples{\dontrun{
data(india)
attach(india)

faccsex <- factor(csex)
facctwin <- factor(ctwin)
faccbirthorder <- factor(cbirthorder)
facmunemployed <- factor(munemployed)
facmreligion <- factor(mreligion)
faccar <- factor(car)

## Estimate a semiparametric additive model averaged model

model <- lm.ma(cheight ~ faccsex + facctwin + faccbirthorder + facmunemployed
                         + facmreligion + faccar + cage + mbmi + medu,
               basis="additive",
               vc=FALSE)

summary(model)

plot(model,plot.data=TRUE)
plot(model,plot.deriv=TRUE)
}}
\keyword{datasets}
