\name{sort.norMix}
\alias{sort.norMix}
\title{Sort Method for "norMix" Objects}
\description{
  Sorting a \code{"norMix"} object (see \code{\link{norMix}}),
  sorts along the \code{mu} values; i.e., for the
  default \code{decreasing = FALSE} the resulting \code{x[,"mu"]} are
  sorted from left to right.
}
\usage{
\method{sort}{norMix}(x, decreasing = FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"norMix"}.}
  \item{decreasing}{logicial indicating if sorting should be up or down.}
  \item{\dots}{further arguments passed to \code{\link{sort}(x[,"mu"],*)}.}
}
\value{
  a \code{"norMix"} object like \code{x}.
}
\examples{
sort(MW.nm9)
stopifnot(identical(MW.nm2, sort(MW.nm2)))
}
\keyword{utilities}
