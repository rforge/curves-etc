\name{qnorMix}
\alias{qnorMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qnorMix(obj, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ ~~Describe \code{obj} here~~ }
  \item{p}{ ~~Describe \code{p} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(obj, p)
{
  if (!is.norMix(obj))
     stop("'obj' must be a 'Normal Mixture' object!")
   mu <- obj[, "mu"]
   sd <- sqrt(obj[, "sig2"])
   n <- nrow(obj)

  ## vectorize in `p' :
  sapply(p, function(p) {
      uniroot(function(l) pnorMix(obj,l) - p,
              interval = range(qnorm(p,mu,sd)))$root
  })
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
