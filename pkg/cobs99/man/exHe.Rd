\name{exHe}
\alias{exHe}
\title{Small Dataset Example of He}
\usage{data(exHe)}
\description{
  The \code{exHe} data frame has 10 rows and 2 columns.  It is an
  example for which \code{\link{smooth.spline}} cannot be used.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{only values 0, 1, and 2.}
    \item{y}{10 randomly generated values}
  }
}
\source{
  Originally found at the bottom of
  \verb{http://ux6.cso.uiuc.edu/~x-he/ftp.html}, the web resource
  directory of Xuming He at the time, say 2006.
}
\examples{
data(exHe)
plot(exHe, main = "He's 10 point example and cobs() fits")
tm <- tapply(exHe$y, exHe$x, mean)
lines(unique(exHe$x), tm, lty = 2)
## 4 warnings :
cH  <- cobsOld(exHe$x, exHe$y, constraint = "increase")
lines(cH$z, cH$fit, col = 3)
cHn <- cobsOld(exHe$x, exHe$y, constraint = "none")
lines(cHn$z, cHn$fit, col = 4)
cHd <- cobsOld(exHe$x, exHe$y, constraint = "decrease")
lines(cHd$z, cHd$fit, col = 5)
}
\keyword{datasets}
