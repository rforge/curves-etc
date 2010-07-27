\name{dnorMix}
\alias{dnorMix}
\alias{dnorMixL}
\alias{dpnorMix}
\title{Normal Mixture Density}
\description{
  Evaluate the density function of the normal mixture specified as
  \code{\link{norMix}} object.
}
\usage{
dnorMix (x, obj,        log = FALSE)

dnorMixL(obj, x = NULL, log = FALSE, xlim = NULL, n = 511)
dpnorMix(x, obj, lower.tail = TRUE)
}
\arguments{
  \item{obj}{an object of class \code{norMix}.}
  \item{x}{numeric vector with abscissa values where to evaluate the
    density (and probability, for \code{dpnorMix()}).  For \code{dnorMixL()}
    by default, when \code{NULL}, it is constructed
    from \code{n} (and \code{xlim} if that is specified).}
  \item{log}{logical indicating \emph{log}-density values should be returned.}
  \item{xlim}{range of abscissa values, used if \code{x == NULL}.  By default,
    xlim is taken as mean plus/minus 3 standard deviations of the normal
    mixture.}
  \item{n}{number of abscissa values to generate if \code{x} is not specified.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dnorMix(x)} returns the numeric vector of density values
  \eqn{f(x)}, logged if \code{log} is TRUE.

  \code{dnorMixL()} returns a list with components
  \item{x}{the abscissa values.}
  \item{y}{the density values \eqn{f(x)} as for \code{dnorMix()}.}

  \code{dpnorMix()} returns a list with components
  \item{d}{the density values \eqn{f(x)} as for \code{dnorMix()}.}
  \item{p}{the probability values \eqn{F(x)} as for \code{pnorMix()}.}
}
\seealso{\code{\link{rnorMix}} for random number generation, and
  \code{\link{norMix}} for the construction and further methods,
  particularly \code{\link{plot.norMix}} which makes use \code{dnorMix}.
}
\examples{
 ff <- dnorMixL(MW.nm7)
 str(ff)
 plot(ff, type = "h", ylim = c(0,1)) # rather use plot(ff, ...)

 x <- seq(-4,5, length=501)
 dp <- dpnorMix(x, MW.nm7)
 lines(x, dp$d, col = "tomato", lwd=3)
 lines(x, dp$p, col = 3, lwd=2)# does not fit y-wise
 stopifnot(all.equal(dp$d, dnorMix(x, MW.nm7), tol=1e-12),
           all.equal(dp$p, pnorMix(x, MW.nm7), tol=1e-12))
}
\keyword{distribution}
