\name{pnorMix}
\encoding{latin1}% Danish oe below
\alias{pnorMix}
\alias{qnorMix}
\title{Normal Mixture Cumulative Distribution and Quantiles}
\description{
  Compute cumulative probabilities or quantiles (the inverse) for a
  normal mixture specified as \code{\link{norMix}} object.
}
\usage{
pnorMix(q, obj, lower.tail = TRUE, log.p = FALSE)

qnorMix(p, obj, lower.tail = TRUE, log.p = FALSE,
        tol = .Machine$double.eps^0.25, maxiter = 1000)
}
\arguments{
  \item{obj}{an object of class \code{norMix}.}
  \item{p}{numeric vector of probabilities.}
  \item{q}{numeric vector of quantiles}.
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{tol, maxiter}{determine the root finding algorithm, see
    \code{\link{uniroot}}.}
}
\details{
  Whereas the distribution function \code{pnorMix} is the trivial sum of
  weighted normal probabilities (\code{\link{pnorm}}), its inverse is
  currently computed numerically using \code{\link{uniroot}} to find
  \code{q} such that \code{pnorMix(obj, q) == p}.
}
\value{
  a numeric vector of the same length as \code{p} or \code{q}, respectively.
}
\author{First version by Erik Jørgensen
  \email{Erik.Jorgensen@agrsci.dk}.}
\seealso{\code{\link{dnorMix}} for the density function.}
\examples{
MW.nm3 # the "strange skew" one
plot(MW.nm3)
## now the cumlative :
x <- seq(-4,4, length=1001)
plot(x, pnorMix(x, MW.nm3), type="l", col=2)
## and some of its inverse :
pp <- seq(.1, .9, by=.1)
plot(qnorMix(pp, MW.nm3), pp)

## The "true" median of a normal mixture:
median.norMix <- function(x) qnorMix(1/2, x)
median.norMix(MW.nm3) ## -2.32
}
\keyword{distribution}
