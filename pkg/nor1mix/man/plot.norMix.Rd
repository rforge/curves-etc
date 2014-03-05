\name{plot.norMix}
\alias{plot.norMix}
\alias{lines.norMix}
\title{Plotting Methods for 'norMix' Objects}
\description{
  The \code{plot} and \code{lines} methods for \code{\link{norMix}}
  objects draw the normal mixture density, optionally additonally with a
  fitted normal density.
}
\usage{
\method{plot}{norMix}(x, type = "l", n = 511, xout = NULL, xlim = NULL,
      xlab = "x", ylab = "f(x)", main = attr(x, "name"), lwd = 1.4,
      p.norm = !p.comp, p.h0 = TRUE, p.comp = FALSE,
      parNorm = list(col = 2, lty = 2, lwd = 0.4),
      parH0 = list(col = 3, lty = 3, lwd = 0.4),
      parComp = list(col= "blue3", lty = 3, lwd = 0.4), \dots)

\method{lines}{norMix}(x, type = "l", n = 511, xout = NULL,
      lwd = 1.4, p.norm = FALSE, parNorm = list(col = 2, lty = 2, lwd = 0.4),
      \dots)
}
\arguments{
  \item{x}{object of class \code{norMix}.}
  \item{type}{character denoting type of plot, see, e.g. \code{\link{lines}}.}
  \item{n}{number of points to generate if \code{xout} is unspecified.}
  \item{xout}{numeric or \code{NULL} giving the abscissae at which to
    draw the density.}
  \item{xlim}{range of x values to use; particularly important if
    \code{xout} is not specified where \code{xlim} is passed to
    \code{\link{dnorMix}} and gets a smart default if unspecified.}
  \item{xlab,ylab}{labels for the x and y axis with defaults.}
  \item{main}{main title of plot, defaulting to the \code{\link{norMix}}
    name.}
  \item{lwd}{line width for plotting with a non-standard default.}
  \item{p.norm}{logical indicating if the normal density with the same
    mean and variance should be drawn as well.}
  \item{p.h0}{logical indicating if the line \eqn{y = 0} should be drawn.}
  \item{p.comp}{logical indicating if the Gaussian components should
    also be drawn individually.}
  \item{parNorm}{graphical parameters for drawing the normal density if
    \code{p.norm} is true.}
  \item{parH0}{graphical parameters for drawing the line \eqn{y = 0} if
    \code{p.h0} is true.}
  \item{parComp}{graphical parameters for drawing the single components if
    \code{p.comp} is true.}
  \item{\dots}{further arguments passed to and from methods.}
}
\author{Martin Maechler}
\seealso{\code{\link{norMix}} for the construction and further methods,
  particularly \code{\link{dnorMix}} which is used here.
}
\examples{
plot(norMix(m=c(0,3), sigma = c(2,1))) # -> var = c(2^2, 1) = c(4, 1)

plot(MW.nm4, p.norm=FALSE, p.comp = TRUE)
## Further examples in  ?norMix and  ?rnorMix
}
\keyword{distribution}
\keyword{hplot}
