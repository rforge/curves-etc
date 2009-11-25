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
        tol = .Machine$double.eps^0.25, maxiter = 1000, traceRootsearch = 0,
        method = c("interpQspline", "interpspline", "eachRoot", "root2"),
        l.interp = 20)
}
\arguments{
  \item{obj}{an object of class \code{norMix}.}
  \item{p}{numeric vector of probabilities.  Note that for all
    \code{method}s but \code{"eachRoot"}, \code{qnorMix(p, *)} works
    with the full vector \code{p}, typically using (inverse)
    interpolation approaches; consequently the result is very
    slightly dependent on \code{p} as a whole.}
  \item{q}{numeric vector of quantiles}.
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{tol, maxiter}{tolerance and maximal number of iterations for the
    root search algorithm, see \code{method} below and \code{\link{uniroot}}.}
  \item{traceRootsearch}{logical or integer in \eqn{\{0,1,2,3\}},
    determining the amount of information printed during root search.}
  \item{method}{a string specifying which algorithm is used for the
        \dQuote{root search}. .......}% <<<< FIXME
  \item{l.interp}{positive integer for \code{method = "interQpspline"}
    or \code{"interpspline"}, determining the number of values in each
    \dQuote{mu-interval}.}
}
\details{
  Whereas the distribution function \code{pnorMix} is the trivial sum of
  weighted normal probabilities (\code{\link{pnorm}}), its inverse,
  \code{qnorMix} is computed numerically: For each \code{p} we search for
  \code{q} such that \code{pnorMix(obj, q) == p}, i.e., \eqn{f(q) = 0}
  for \code{f(q) := pnorMix(obj, q) - p}.  This is a root finding
  problem which can be solved by \code{\link{uniroot}(f, lower,upper,*)}.
  If \code{length(p) <= 2} or \code{method = "eachRoot"}, this happens
  one for one for the \emph{sorted} p's.  Otherwise, we start by doing
  this for the outermost non-trivial (\eqn{0 < p < 1}) values of p.

  %% "root2" may be abandoned, so I'm too lazy to describe it now --- FIXME
  %%
  For \code{method = "interQpspline"} or \code{"interpspline"}, we now compute
  \code{p. <- \link{pnorMix}(q., obj)} for values \code{q.} which are a grid
  of length \code{l.interp} in each interval \eqn{[q_j,q_{j+1}]} and
  then use \emph{montone} inverse interpolation
  (\code{\link{splinefun}(q., p., method="monoH.FC")}) plus
  a few (maximally \code{maxiter}, typically one!) Newton steps.
  The default, \code{"interQpspline"}, additionally logit-transforms the
  \code{p.} values to make the interpolation more linear.
  This method is faster, particularly for large \code{length(p)}.
}
\value{
  a numeric vector of the same length as \code{p} or \code{q}, respectively.
}
\author{Very first version (for length-1 \code{p,q}) by
  Erik Jørgensen \email{Erik.Jorgensen@agrsci.dk}.
}
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
