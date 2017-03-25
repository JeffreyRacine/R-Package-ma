\name{lm.ma}
\alias{lm.ma}
\title{
Model Averaging
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lm.ma(y = NULL, X = NULL, X.eval = NULL, basis = c("glp", "tensor", "additive"), compute.deriv = TRUE, p.max = 25, method = c("jma", "mma"), ma.weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{X.eval}{
%%     ~~Describe \code{X.eval} here~~
}
  \item{basis}{
%%     ~~Describe \code{basis} here~~
}
  \item{compute.deriv}{
%%     ~~Describe \code{compute.deriv} here~~
}
  \item{p.max}{
%%     ~~Describe \code{p.max} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{ma.weights}{
%%     ~~Describe \code{ma.weights} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
X.eval <- expand.grid(x.seq,x.seq)

out <- lm.ma(y~x1+x2)

z <- matrix(out$fitted,n.eval,n.eval)

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
