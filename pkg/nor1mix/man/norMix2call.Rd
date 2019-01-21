\name{norMix2call}
\alias{norMix2call}
\alias{as.expression.norMix}
\alias{as.function.norMix}
\title{Transform "norMix" object into Call, Expression or Function}
\description{
  E.g., for taking symbolic derivatives, it may be useful to get an \R
  \code{\link{call}}, \code{\link{expression}}, or \code{\link{function}}
  in / of \code{x} from a specific \code{"\link{norMix}"} object.
}
\usage{
norMix2call(obj, oneArg = TRUE)
\method{as.expression}{norMix}(x, oneArg = TRUE, \dots)
\method{as.function}{norMix}  (x, oneArg = TRUE, envir = parent.frame(), \dots)
}
\arguments{
  \item{obj, x}{an \R object of \code{\link{class}} \code{"\link{norMix}"}.}
  \item{oneArg}{\code{\link{logical}} specifying if expressions of the form
    \code{dnorm((x - mu)/sig)} should be used, i.e. \bold{one Arg}ument
    only, instead of \code{dnorm(x, mu, sig)}.}
  \item{envir}{an \code{\link{environment}}; often the default is perfect.}
  \item{\dots}{potentially further arguments (not used in any examples yet).}
}
\value{
  according to the function used, an \R  \sQuote{language} object, i.e., a
  \code{\link{call}}, \code{\link{expression}}, or \code{\link{function}}, respectively.
}
\author{Martin Maechler}
\seealso{
  \code{\link{norMix}}.  Note that \code{\link{deriv}()} currently only
  works correctly in case of the default \code{oneArg = TRUE}.
}
\examples{
(cMW2 <- norMix2call(MW.nm2))
deriv(cMW2, "x")

(fMW1 <- as.function  (MW.nm1))
(eMW3 <- as.expression(MW.nm3))
stopifnot(is.call      (cMW2),  is.call(norMix2call(MW.nm2, FALSE)),
          is.function  (fMW1),  is.function  (as.function  (MW.nm4)),
          is.expression(eMW3),  is.expression(as.expression(MW.nm5))
)
}
\keyword{math}
