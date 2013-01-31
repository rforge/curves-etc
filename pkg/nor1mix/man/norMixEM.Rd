\name{norMixEM}
\alias{norMixEM}
\title{
  EM Estimation of Mixtures of Univariate Normal Distributions
}
\description{
  This function estimates the parameters of a univariate normal mixture
  using the EM algorithm.
}
\usage{
norMixEM(x, m, name=NULL, name=NULL, iter.max = 100,
         tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{x}{numeric: the data for which the parameters are to be estimated.}
  \item{m}{integer or factor: If \code{m} has length 1 it specifies the
    number of mixture components, otherwise it is taken to be a vector
    of initial cluster assignments, see details below.}
  \item{name}{character: passed to \code{\link{norMix}}.}
  \item{iter.max}{integer: maximum number of EM iterations.}
  \item{tol}{numeric: EM iterations stop if relative changes of the
    log-likelihood are smaller than \code{tol}.}
}
\details{
  EM estimation of univariate mixtures can be very sensitive to
  initialization. By default \code{norMixEM} \code{\link{cut}}s the data
  into groups of approximately equal size. See examples below for other
  initialization possibilities.
}
\value{
  An object of class \code{\link{norMix}}.
}
\author{
  Friedrich Leisch
}
\examples{
ex <- norMix(mu = c(-1,2,5), sig2=c(1, 0.5, 3))
plot(ex)

x <- rnorMix(100, ex)
plot(density(x))
rug(x)

## EM estimation may fail depending on random sample
ex1 <- norMixEM(x, 3)
ex1
plot(ex1)

## initialization by cut() into intervals of equal length:
ex2 <- norMixEM(x, cut(x, 3))
ex2

## initialization by kmeans():
ex3 <- norMixEM(x, kmeans(x, 3)$cluster)
ex3

## real data
data(faithful)
plot(density(faithful$waiting))
ex4 <- norMixEM(faithful$waiting, 2)
ex4
plot(ex4)
}
