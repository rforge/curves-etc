\name{varNPreg}
\alias{varNPreg}
\title{Nonparametric Variance Estimator}
\description{
  Estimates the error variance \eqn{\sigma^2} nonparametrically in the model
  \deqn{Y_i = m(x_i) + E_i,} where
  \eqn{E_i \sim (0,\sigma^2)}{E_i ~ (0,sigma^2)}, i.i.d.

  Computes leave-one-out residuals (local linear approximation followed by
  reweighting) and their variance.
}
\usage{
varNPreg(x, y)
}
\arguments{
  \item{x}{abscissae values, ordered increasingly.}
  \item{y}{observations at \code{y[i]} at \code{x[i]}.}
}
\note{
  This is an \R interface to the \code{resest} Fortran subroutine, used
  in \code{\link{lokerns}} and \code{\link{glkerns}}, see their help
  pages for references and context.
}
\value{
  A list with components
  \item{res}{numeric; residuals at \code{x[]} of length \code{n}.}
  \item{snr}{explained variance of the true curve}
  \item{sigma2}{estimation of residual variance, \eqn{\hat\sigma^2}{sigma^2}.}
}
\author{Martin Maechler}
\seealso{\code{\link{lokerns}}, \code{\link{glkerns}}.}
\examples{
 x <- sort(runif(100))
 y <- sin(pi*x) + rnorm(100)/10
 str(ve <- varNPreg(x,y))
}
\keyword{nonparametric}
\keyword{regression}
