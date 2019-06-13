\name{llnorMix}
\title{Likelihood, Parametrization and EM-Steps For 1D Normal Mixtures}
\alias{llnorMix}
\alias{par2norMix}
\alias{nM2par}
\alias{estep.nm}
\alias{mstep.nm}
\alias{emstep.nm}
\description{
  These functions work with an almost unconstrained parametrization of
  univariate normal mixtures.
  \describe{
    \item{\code{llnorMix(p, *)}}{computes the log likelihood,}
    \item{\code{obj <- par2norMix(p)}}{maps parameter vector \code{p} to
      a \code{\link{norMix}} object \code{obj},}
    \item{\code{p <- nM2par(obj)}}{maps from a \code{\link{norMix}}
      object \code{obj} to parameter vector \code{p},}
  }
  where \code{p} is always a parameter vector in our parametrization.

  Partly for didactical reasons, the following functions provide the
  basic ingredients for the EM algorithm (see also
  \code{\link{norMixEM}}) to parameter estimation:
  \describe{
    \item{\code{estep.nm(x, obj, p)}}{computes 1 E-step for the data
      \code{x}, given \emph{either} a \code{"norMix"} object \code{obj}
      or parameter vector \code{p}.}

    \item{\code{mstep.nm(x, z)}}{computes 1 M-step for the data
      \code{x} and the probability matrix \code{z}.}

    \item{\code{emstep.nm(x, obj)}}{computes 1 E- and 1 M-step for the data
      \code{x} and the \code{"norMix"} object \code{obj}.}
  }
  where again, \code{p} is a parameter vector in our parametrization,
  \code{x} is the (univariate) data, and \code{z} is a \eqn{n \times
    m}{n * m} \code{\link{matrix}} of (posterior) conditional
  probabilities, and \eqn{\theta} is the full parameter vector of the
  mixture model.
}
\usage{
llnorMix(p, x, m = (length(p) + 1)/3, trafo = c("clr1", "logit"))

par2norMix(p, trafo = c("clr1", "logit"), name = )

nM2par(obj, trafo = c("clr1", "logit"))

 estep.nm(x, obj, par)
 mstep.nm(x, z)
emstep.nm(x, obj)
}
\arguments{
  \item{p, par}{numeric vector: our parametrization of a univariate normal
    mixture, see details.}
  \item{x}{numeric: the data for which the likelihood is to be computed.}
  \item{m}{integer number of mixture components; this is not to be
    changed for a given \code{p}.}
  \item{trafo}{\code{\link{character}} string specifying the transformation
    of the component weight \code{w} \eqn{m}-vector (mathematical notation
    in \code{\link{norMix}}: \eqn{\pi_j, j=1,\dots,m}) to an
    \eqn{m-1}-dimensional unconstrained parameter vector in our
    parametrization.  \code{"logit"} has been hard-wired upto \pkg{nor1mix}
    version 1.2-3, and has been replaced \emph{as default} in 2019 for \pkg{nor1mix}
    version 1.2-4 by \code{"clr1"} which is more symmetric and basically
    Aitchinson's \bold{c}entered \bold{l}og \bold{r}atio, see also CRAN package
    \CRANpkg{compositions}' \code{clr()}.}
  \item{name}{(for \code{par2norMix()}:) a name for the \code{"norMix"}
    object that is returned, uses a smart default.}
  \item{obj}{a \code{"norMix"} object, see \code{\link{norMix}}.}
  \item{z}{a \eqn{n \times m}{n * m} \code{\link{matrix}} of (posterior)
    conditional probabilities,
    \eqn{z_{ij}= P(x_i \in C_j | \theta)}{z[i,j]= P(x[i] in C_j | theta)},
    where \eqn{C_j} denotes the \eqn{j}-th group (\dQuote{cluster}).}
}
\details{
  We use a parametrization of a (finite) \eqn{m}-component univariate normal mixture which
  is particularly apt for likelihood maximization, namely, one whose
  parameter space is \emph{almost} a full \eqn{\mathbf{I\hskip-0.22em R}^M}{R^M},
  \eqn{M = 3m-1}.

  For an \eqn{m}-component mixture,
  we map to and from a parameter vector \eqn{\theta} (\code{== p} as \R-vector)
  of length \eqn{3m-1}.  For mixture density
  \deqn{\sum_{j=1}^m \pi_j \phi((t - \mu_j)/\sigma_j)/\sigma_j,}{
    sum[j=1..m]  pi[j] phi((t - mu[j])/s[j])/s[j],}
  we transform the \eqn{\pi_j}{pi[j]} (for \eqn{j \in 1,\dots,m}{j =
  1,\dots,m}) via the transform specified by \code{trafo} (see below),
  and log-transform the \eqn{\sigma_j}{s[j]}.  Consequently, \eqn{\theta} is
  partitioned into
  \describe{
    \item{\code{p[ 1:(m-1)]}: }{For
      \describe{
	\item{\code{trafo = "logit"}:}{
	  \code{p[j]}\eqn{=
	    \mbox{logit}(\pi_{j+1})}{= logit(pi[j+1])} and
	  \eqn{\pi_1}{pi[1]} is given implicitly as
	  \eqn{\pi_1 = 1-\sum_{j=2}^m \pi_j}{pi[1] = 1 - sum[j=2..m] pi[j]}.}

	\item{\code{trafo = "clr1"}:}{(\bold{c}entered \bold{l}og
	  \bold{r}atio, omitting 1st element):
	  Set \eqn{\ell_j := \ln(\pi_j)}{l_j := ln(pi[j])} for \eqn{j = 1, \dots, m}, and
	  \code{p[j]}\eqn{= \ell_{j+1} - 1/m\sum_{j'=1}^m \ell_{j'}}{= l_{j+1} - mean_j' l_j'}
	  for \eqn{j = 1, \dots, m-1}.}
      }
    }

    \item{\code{p[ m:(2m-1)]}: }{\code{p[m-1+ j]}\eqn{= \mu_j}, for j=1:m.}

    \item{\code{p[2m:(3m-1)]}: }{\code{p[2*m-1+ j]}%
      \eqn{= \log(\sigma_j)}{= log(s[j])}, i.e.,
      \eqn{\sigma_j^2 = exp(2*p[.+j])}.}
  }
}
\value{
  \code{llnorMix()} returns a number, namely the log-likelihood.

  \code{par2norMix()} returns \code{"norMix"} object, see \code{\link{norMix}}.

  \code{nM2par()} returns the parameter vector \eqn{\theta} of length \eqn{3m-1}.

  \code{estep.nm()} returns \code{z}, the matrix of (conditional) probabilities.

  \code{mstep.nm()} returns the model parameters as a
  \code{\link{list}} with components \code{w}, \code{mu}, and
  \code{sigma}, corresponding to the arguments of \code{\link{norMix}()}.
  (and see the 'Examples' on using \code{\link{do.call}(norMix, *)} with it.)

  \code{emstep.nm()} returns an \emph{updated} \code{"norMix"} object.
}
%%\references{ MM should write a report ...}
\author{Martin Maechler}
\seealso{\code{\link{norMix}}, \code{\link{logLik}}.
  Note that the log likelihood of a \code{"norMix"} object
  is directly given by \code{sum(\link{dnorMix}(x, obj, log=TRUE))}.

  To fit, using the EM algorithm, rather use \code{\link{norMixEM}()}
  than the \code{e.step}, \code{m.step}, or \code{em.step} functions.

  Note that direct likelihood maximization, i.e., MLE, is typically
  considerably more efficient than the EM, and typically converges well
  with our parametrization, see \code{\link{norMixMLE}}.
}
\examples{
(obj <- MW.nm10) # "the Claw" -- m = 6 components
length(pp <- nM2par(obj)) # 17 == (3*6) - 1
par2norMix(pp)
## really the same as the initial 'obj' above

## Log likelihood (of very artificial data):
llnorMix(pp, x = seq(-2, 2, length=1000))
set.seed(47)## of more realistic data:
x <- rnorMix(1000, obj)
llnorMix(pp, x)

## Consistency check :  nM2par()  and  par2norMix()  are inverses
all.EQ <- function(x,y, tol = 1e-15, ...) all.equal(x,y, tolerance=tol, ...)
stopifnot(all.EQ(pp, nM2par(par2norMix(pp))),
          all.EQ(obj, par2norMix(nM2par(obj)),
                    check.attributes=FALSE),
          ## Direct computation of log-likelihood:
          all.EQ(sum(dnorMix(x, obj, log=TRUE)),
                    llnorMix(pp, x))  )

## E- and M- steps : ------------------------------
rE1 <- estep.nm(x, obj)
rE2 <- estep.nm(x, par=pp) # the same as rE1
z <- rE1
str( rM <-  mstep.nm(x, z))
   (rEM <- emstep.nm(x, obj))

stopifnot(all.EQ(rE1, rE2),
          all.EQ(rEM, do.call(norMix, c(rM, name=""))))
}
\keyword{models}

