\name{norMixFit}
\alias{norMixEM}
\alias{norMixMLE}
\title{EM and MLE Estimation of Univariate Normal Mixtures}
\description{
  These functions estimate the parameters of a univariate (finite) normal
  mixture using the EM algorithm or Likelihood Maximimization via
  \code{\link{optim}(.., method = "BFGS")}.
}
\usage{
norMixEM(x, m, name = NULL, sd.min =  1e-07* diff(range(x))/m,
         maxiter = 100, tol = sqrt(.Machine$double.eps), trace = 1)

norMixMLE(x, m, name = NULL, %% sd.min = 1e-7 * diff(range(x))/m,
          maxiter = 100, tol = sqrt(.Machine$double.eps), trace = 2)
}
\arguments{
  \item{x}{numeric: the data for which the parameters are to be estimated.}
  \item{m}{integer or factor: If \code{m} has length 1 it specifies the
    number of mixture components, otherwise it is taken to be a vector
    of initial cluster assignments, see details below.}
  \item{name}{character, passed to \code{\link{norMix}}.  The default,
    \code{NULL}, uses \code{\link{match.call}()}.}
  \item{sd.min}{number: the minimal value that the normal components'
    standard deviations (\code{sd}) are allowed to take.  A warning is
    printed if some of the final \code{sd}'s are this boundary.}
  \item{maxiter}{integer: maximum number of EM iterations.}
  \item{tol}{numeric: EM iterations stop if relative changes of the
    log-likelihood are smaller than \code{tol}.}
  \item{trace}{integer (or logical) specifying if the iterations should be
    traced and how much output should be produced.  The default,
    \code{1} prints a final one line summary, where \code{trace = 2}
    produces one line of output per iteration.}
}
\details{
  Estimation of univariate mixtures can be very sensitive to
  initialization.  By default, \code{norMixEM} and \code{norMixLME}
  \code{\link{cut}} the data into \eqn{m} groups of approximately equal
  size.  See examples below for other initialization possibilities.

  The EM algorithm consists in repeated application of E- and M- steps
  until convergence.  Mainly for didactical reasons, we also provide the
  functions \code{\link{estep.nm}}, \code{mstep.nm}, and
  \code{\link{emstep.nm}}.

  The MLE, Maximum Likelihood Estimator, maximizes the likelihood using
  \code{\link{optim}}, using the same advantageous parametrization as
  \code{\link{llnorMix}}.
}
\value{
  An object of class \code{\link{norMix}}.
}
\author{
  EM: Friedrich Leisch, originally; Martin Maechler vectorized it in
  \eqn{m}, added \code{trace} etc.

  MLE: M.Maechler
}
\examples{
ex <- norMix(mu = c(-1,2,5), sig2 = c(1, 0.5, 3))
plot(ex, col="gray", p.norm=FALSE)

x <- rnorMix(100, ex)
lines(density(x))
rug(x)

## EM estimation may fail depending on random sample
ex1 <- norMixEM(x, 3, trace=2) #-> warning (sometimes)
ex1
plot(ex1)

## initialization by cut() into intervals of equal length:
ex2 <- norMixEM(x, cut(x, 3))
ex2

## initialization by kmeans():
k3 <- kmeans(x, 3)$cluster
ex3 <- norMixEM(x, k3)
ex3

## Now, MLE instead of EM:
exM <- norMixMLE(x, k3, tol = 1e-12, trace=4)
exM

## real data
data(faithful)
plot(density(faithful$waiting, bw = "SJ"), ylim=c(0,0.044))
rug(faithful$waiting)

(nmF <- norMixEM(faithful$waiting, 2))
lines(nmF, col=2)
## are three components better?
nmF3 <- norMixEM(faithful$waiting, 3, maxiter = 200)
lines(nmF3, col="forestgreen")
}
