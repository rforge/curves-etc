\name{llnorMix}
\title{Likelihood and Parametrization of 1D Normal Mixtures}
\alias{llnorMix}
\alias{par2norMix}
\alias{nM2par}
\description{
  These functions work with an \emph{unconstrained} parametrization of
  univariate normal mixtures.
}
\usage{
llnorMix(p, x, m = (length(p) + 1)/3)

par2norMix(p, name = sprintf("{from \%s}", deparse(substitute(p))))
nM2par(obj)
}
\arguments{
  \item{p}{numeric vector: our parametrization of a univariate normal
    mixture, see details.}
  \item{x}{numeric: the data for which the likelihood is to be computed.}
  \item{m}{integer number of mixture components; this is not to be
    changed for a given \code{p}.}
  \item{name}{(for \code{par2norMix()}:) a name for the \code{"norMix"}
    object that is returned.}
  \item{obj}{a \code{"norMix"} object, see \code{\link{norMix}}.}
}
\details{
  We use a parametrization of a (finite) univariate normal mixture which
  is particularly apt for likelihood maximization, namely, one whose
  parameter space is a full \eqn{\mathbf{I\hskip-0.22em R}^m}{R^m},
  \eqn{m = 3k-1}.

  For a \eqn{k}-component mixture,
  we map (1:1) a parameter vector \eqn{\theta} (\R-vector \code{p}) of
  length \eqn{3k-1}. For mixture density
  \deqn{\sum_{j=1}^k \pi_j \phi((t - \mu_j)/\sigma_j),}{
    sum[j=1..k]  pi[j] phi((t - mu[j])/s[j]),}
  we logit-transform the \eqn{\pi_j}{pi[j]} (for \eqn{j \ge 2}{j >= 2})
  and log-transform the \eqn{\sigma_j}{s[j]}, such that \eqn{\theta} is
  partitioned into
  \describe{
    \item{\code{p[ 1:(k-1)]}: }{\code{p[j]}\eqn{=
	\mbox{logit}(\pi_{j+1})}{= logit(pi[j+1])} and
      \eqn{\pi_1}{pi[1]} is given implicitly as
      \eqn{\pi_1 = 1-\sum_{j=2}^k \pi_j}{pi[1] = 1 - sum[j=2..k] pi[j]}.}

    \item{\code{p[ k:(2k-1)]}: }{\code{p[k-1+ j]}\eqn{= \mu_j}, for j=1:k.}

    \item{\code{p[2k:(3k-1)]}: }{\code{p[2*k-1+ j]}%
      \eqn{= \log(\sigma_j)}{= log(s[j])}, i.e.,
      \eqn{\sigma_j^2 = exp(2*p[.+j])}.}
  }
}
\value{
  \code{llnorMix()} returns a number, namely the log-likelihood.

  \code{par2norMix()} returns \code{"norMix"} object, see \code{\link{norMix}}.

  \code{nM2par()} returns the parameter vector \eqn{\theta} of length
  \eqn{3k-1}.
}
%%\references{ MM should write a report ...}
\author{Martin Maechler}
\seealso{\code{\link{norMix}}, \code{\link{logLik}}
}
\examples{
(obj <- MW.nm10) # "the Claw" -- m = 6 components
length(pp <- nM2par(obj)) # 17 == (3*6) - 1
par2norMix(pp)

## Log likelihood (of very artificial data):
llnorMix(pp, x = seq(-2, 2, by=.1))

## Consistency check :
stopifnot(all.equal(pp, nM2par(par2norMix(pp)), tol= 1e-15),
          all.equal(obj, par2norMix(nM2par(obj)),
                    check.attributes=FALSE, tol=1e-15))
}
\keyword{models}

