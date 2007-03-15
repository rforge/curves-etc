\name{dnorMix}
\alias{dnorMix}
\alias{dnorMixL}
\title{Normal Mixture Density}
\description{
  Evaluate the density function of the normal mixture specified as
  \code{\link{norMix}} object.
}
\usage{
dnorMix (x, obj,        log = FALSE)
dnorMixL(obj, x = NULL, log = FALSE, xlim = NULL, n = 511)
}
\arguments{
  \item{obj}{an object of class \code{norMix}.}
  \item{x}{numeric vector with abscissa values where to evaluate the
    density.  For \code{dnorMixL()} by default, when \code{NULL}, it is
    constructed from \code{n} (and \code{xlim} if that is specified).}
  \item{log}{logical indicating \emph{log}-density values should be returned.}
  \item{xlim}{range of abscissa values, used if \code{x == NULL}.  By default,
    xlim is taken as mean plus/minus 3 standard deviations of the normal
    mixture.}
  \item{n}{number of abscissa values to generate if \code{x} is not specified.}
}
\value{
  \code{dnorMix(x)} returns the numeric vector of density values
  \eqn{f(x)}, logged if \code{log} is TRUE.

  \code{dnorMixL()} returns a list with components
  \item{x}{the abscissa values.}
  \item{y}{the density values \eqn{f(x)} as for \code{dnorMix()}.}
}
\seealso{\code{\link{rnorMix}} for random number generation, and
  \code{\link{norMix}} for the construction and further methods,
  particularly \code{\link{plot.norMix}} which makes use \code{dnorMix}.
}
\examples{
 ff <- dnorMixL(MW.nm7)
 str(ff)
 plot(ff, type = "h") # rather use plot(ff, ...)
}
\keyword{distribution}
