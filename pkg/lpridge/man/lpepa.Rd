\name{lpepa}
\encoding{latin1}% for a"
\alias{lpepa}
\title{Local polynomial regression fitting with Epanechnikov weights}
\description{
  Fast and stable algorithm for nonparametric estimation of regression
  functions and their derivatives via \bold{l}ocal \bold{p}olynomials with
  \bold{Epa}nechnikov weight function.
}
\usage{
lpepa(x, y, bandwidth, deriv = 0, n.out = 200, x.out = NULL,
      order = deriv+1, mnew = 100, var = FALSE)
}
\arguments{
  \item{x}{vector of design points, not necessarily ordered.}
  \item{y}{vector of observations of the same length as \code{x}.}
  \item{bandwidth}{bandwidth(s) for nonparametric estimation.  Either a
    number or a vector of the same length as \code{x.out}.}
  \item{deriv}{order of derivative of the regression function to be
    estimated; defaults to \code{deriv = 0}.}
  \item{n.out}{number of output design points where the function has to
    be estimated.  The default is \code{n.out=200}.}
  \item{x.out}{vector of output design points where the function has to
    be estimated.  The default value is an equidistant grid of \code{n.out}
    points from min(x) to max(x).}
  \item{order}{integer, order of the polynomial used for local
    polynomials.  Must be \eqn{\le 10}{<= 10} and defaults to
    \code{order = deriv+1}.}
  \item{mnew}{integer forcing to restart the algorithm after \code{mnew}
    updating steps.
    The default is \code{mnew = 100}.  For \code{mnew = 1} you get a numerically
    \dQuote{super-stable} algorithm (see reference SBE\&G below).}
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
\details{More details are described in the first reference SBE\&G (1994)
  below.  In S\&G, a bad finite sample behaviour of local polynomials
  for random designs was found.
  For practical use, we therefore propose local polynomial regression
  fitting with ridging, as implemented in the function
  \code{\link{lpridge}}.  In \code{lpepa}, several parameters described
  in SBE\&G are fixed either in the fortran
  routine or in the \R-function.  There, you find comments how to change them.

  %% MM: FIXME -- make this more transparent and user configurable

  For \code{var=TRUE}, the variance of the estimator proportional to the
  residual variance is computed, i.e., the exact finite sample variance of the
  regression estimator is \code{var(est) = est.var * sigma^2}.
}
\references{
  See also \url{http://www.unizh.ch/biostat/} under
  \file{Manuscripts} etc.

  - Numerical stability and computational speed:
  \cr
  B. Seifert, M. Brockmann, J. Engel and T. Gasser (1994)
  Fast algorithms for nonparametric curve estimation.
  \emph{J. Computational and Graphical Statistics} \bold{3}, 192--213.

  - Statistical properties:
  \cr
  Seifert, B. and Gasser, T. (1996)
  Finite sample variance of local polynomials: Analysis and solutions.
  \emph{J. American Statistical Association} \bold{91}(433), 267--275.

  Seifert, B. and Gasser, T. (2000)
  Data adaptive ridging in local polynomial
  regression.  \emph{J. Computational and Graphical Statistics} \bold{9},
  338--360.

  Seifert, B. and Gasser, T. (1998)
  Ridging Methods in Local Polynomial Regression.
  in: S. Weisberg (ed), \emph{Dimension Reduction, Computational Complexity,
  and Information}, Vol.\bold{30} of Computing Science \& Statistics,
  Interface Foundation of North America, 467--476.

  Seifert, B. and Gasser, T. (1998)
  Local polynomial smoothing.
  in: \emph{Encyclopedia of Statistical Sciences},
  Update Vol.\bold{2}, Wiley, 367--372.

  Seifert, B., and Gasser, T. (1996)
  Variance properties of local polynomials and ensuing
  modifications. in: \emph{Statistical Theory and Computational Aspects
    of Smoothing}, W. Härdle, M. G. Schimek (eds), Physica, 50--127.
}
\seealso{\code{\link{lpridge}}, and also \code{\link{lowess}} and
  \code{\link{loess}} which do local linear and quadratic regression
  quite a bit differently.
}
\examples{
data(cars)
attach(cars)

epa.sd <- lpepa(speed,dist, bandw=5)		# local polynomials

plot(speed, dist, main = "data(cars) & lp epanechnikov regression")
lines(epa.sd$x.out, epa.sd$est,  col="red")
lines(lowess(speed,dist, f= .5), col="orange")
detach()
}
\keyword{smooth}

