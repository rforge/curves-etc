\name{lpepa}
\alias{lpepa}
\title{Local polynomial regression fitting with Epanechnikov weights}
\usage{
lpepa(x, y, bandwidth, deriv = 0, n.out = 200, x.out = NULL,
      order=NULL, mnew=100, var=F)
}
\arguments{
  \item{x}{vector of design points, not necessarily ordered.}
  \item{y}{vector of observations of the same length as x.}
  \item{bandwidth}{bandwidth for nonparametric estimation.  Either a
    number or a vector of the same length as \code{x.out}.}
  \item{deriv}{order of derivative of the regression function to be
    estimated. The default value is deriv=0.}
  \item{n.out}{number of output design points where the function has to
    be estimated.  The default value is n.out=200.}
  \item{x.out}{vector of output design points where the function has to
    be estimated.  The default value is an equidistant grid of \code{n.out}
    points from min(x) to max(x).}
  \item{order}{order of the polynomial used for local polynomials.  The
    default value is order=deriv+1.}
  \item{mnew}{force to restart the algorithm after mnew updating steps.
    The default value is mnew=100.  For mnew=1 you get a numerically
    "super-stable" algorithm (see reference SBE&G below).}
  \item{var}{logical flag: if \code{TRUE}, the variance of the estimator
    proportional to the residual variance is computed (see details).}
}
\value{a list including used parameters and estimator.
  \item{x}{vector of ordered design points.}
  \item{y}{vector of observations ordered according to x.}
  \item{bandwidth}{vector of bandwidths actually used for nonparametric
    estimation.}
  \item{deriv}{order of derivative of the regression function estimated.}
  \item{x.out}{vector of ordered output design points.}
  \item{order}{order of the polynomial used for local polynomials.}
  \item{mnew}{force to restart the algorithm after mnew updating steps.}
  \item{var}{logical flag: whether the variance of the estimator was computed.}
  \item{est}{estimator of the derivative of order deriv of the
    regression function.}
  \item{est.var}{estimator of the variance of est (proportional to
    residual variance).}
}
\description{
  Fast and stable algorithm for nonparametric estimation of regression
  functions and their derivatives via local polynomials with
  Epanechnikov weight function.}
\details{More details are described in the reference SBE\\&G below.  In
  S\\&G a bad finite
  sample behaviour of local polynomials for random design was found.
  For practical use we propose local polynomial regression fitting with
  ridging, as implemented in the function "lpridge".  In lpepa, several
  parameters described in SBE\\&G are fixed either in the fortran
  routine or in the Splus-function. There, you find comments how to
  change them.
  For var=T, the variance of the estimator proportional to the residual
  variance is computed, i.e. the exact finite sample variance of the
  regression estimator is var(est) = est.var*sigma^2.
}
\references{
  - numerical stability and computational speed:\cr
  B. Seifert, M. Brockmann, J. Engel & T. Gasser (1994).
  Fast algorithms for nonparametric curve estimation.
  J. Computational and Graphical Statistics 3, 192-213.

  - statistical properties:\cr
  B. Seifert & T. Gasser (1994).
  Finite sample variance of local polynomials: Analysis and solutions.
  submitted (also on WWW: http://www.unizh.ch/biostat).
}
\seealso{\code{\link{lpridge}}
\examples{
data(cars)
attach(cars)

epa.sd <- lpepa(speed,dist, bandw=5)		# local polynomials

plot(speed, dist, main = "data(cars) & lp epanechnikov regression")
lines(epa.sd$x.out, epa.sd$est, col="red")
lines(lowess(speed,dist, f= .5),col="orange")
detach()
}
\keyword{smoothing}

