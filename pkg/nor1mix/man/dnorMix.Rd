\name{dnorMix}
\alias{dnorMix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dnorMix(obj, x = NULL, xlim = NULL, n = 511, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{obj}{ ~~Describe \code{obj} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of `comp1'}
  \item{comp2 }{Description of `comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(obj, x = NULL, xlim = NULL, n = 511,...)
{
  ## Purpose: density evaluation for "norMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Arguments: obj: Normal Mixture object;  x: abscissa values where to eval
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  if(is.null(x)) {
    if(is.null(xlim)) { ##-- construct "reasonable" abscissa values
      xlim <- mean.norMix(obj) + c(-3,3)*sqrt(var.norMix(obj))
    }
    x <- seq(xlim[1], xlim[2], length = n)
  }
  m <- m.norMix(obj) #-- number of components
  y <- numeric(length(x))
  w <- obj[,"w"]; mu <- obj[,"mu"]; sig2 <- obj[,"sig2"]
  for(j in 1:m)
    y <- y + w[j] * dnorm(x, mean = mu[j], sd = sqrt(sig2[j]))
  list(x = x, y = y)
  }
}
\keyword{distribution}
