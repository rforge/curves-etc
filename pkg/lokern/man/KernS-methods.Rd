\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\name{KernS-methods}
\title{Methods for ("KernS" classed) Results of lokerns() and glkerns()}
\alias{fitted.KernS}
\alias{lines.KernS}
\alias{plot.KernS}
\alias{predict.KernS}
\alias{print.KernS}
\alias{residuals.KernS}
\description{
  Methods for results of \code{\link{glkerns}()} and
  \code{\link{lokerns}()} which are of (S3) class \code{"KernS"}.
}
\usage{
\S3method{fitted}{KernS}(object, \dots)
\S3method{plot}{KernS}(x, type = "l", lwd = 2.5, col = 3, \dots)
\S3method{predict}{KernS}(object, x, deriv = object[["deriv"]],
        korder = deriv+2, trace.lev = 0, \dots)
\S3method{print}{KernS}(x, digits = getOption("digits"), \dots)
\S3method{residuals}{KernS}(object, \dots)
}
\arguments{
  \item{x,object}{an \R object, of S3 class \code{"KernS"}, typically
    result either from \code{\link{glkerns}()} or \code{\link{lokerns}()}.}
  \item{type, lwd, col}{arguments for \code{plot()} \emph{only} for the
    case when \code{x$deriv} is \emph{not} 0.}
  \item{deriv}{integer, \eqn{\ge 0}{>= 0}, specifiying order of
    derivative that should be predicted.}
  \item{korder}{nonnegative integer giving the kernel order; see
    \code{\link{lokerns}} or \code{\link{glkerns}}.}
  \item{digits}{number of significant digits, see \code{\link{print}}.}
  \item{trace.lev}{integer; level of tracing of Fortran level
    computations; see \code{\link{lokerns}}.}
  \item{\dots}{potentially further arguments passed to and from
    methods.  For the \code{plot(*, deriv=0)} method, these are passed to
    \code{\link[sfsmisc]{plotDS}} from package \CRANpkg{sfsmisc}.}
}
\details{
  Note that \code{fitted()} and \code{residuals()} rely on
  \code{x.inOut} having been true or \code{x.out} having contained the
  data \code{x}, in the \code{\link{lokerns}} or \code{\link{glkerns}}
  call.

  The \code{plot()} method calls \code{\link{plotDS}} from
  package \CRANpkg{sfsmisc}.

  \code{predict(object, x, deriv)} when either some \code{x} are not in
  \code{x.out} or \code{deriv} is not 0, basically recalls the original
  \code{\link{lokerns}} or \code{\link{glkerns}} function (keeping the
  bandwidths for \code{lokerns}).
}
\value{
  (differing, depending on the generic function)%% FIXME
}
\seealso{
 \code{\link{glkerns}}, \code{\link{lokerns}}.
}
\examples{
## "interesting" artificial data:
set.seed(47)
x <- sort(round(10*runif(250),2))
fx <- 5 - x/2 + 3*exp(-(x-5)^2)
y <- fx + rnorm(fx)/4
plot(x,y)
lof <- lokerns(x,y, trace=1)# tracing the phases inside the Fortran code
plot(lof)
plot(lof, cex = 1/4)# maybe preferable
## Simpler, using the lines() method:
plot(x,y); lines(lof, lwd=2, col=2)

qqnorm(residuals(lof)) # hmm... overfitting?
stopifnot(all.equal(y, fitted(lof) + residuals(lof), tolerance = 1e-15),
          predict(lof)$y == fitted(lof))
lof$iter # negative ?
tt <- seq(0, 10, by=1/32)
p0 <- predict(lof, x=tt)
p1 <- predict(lof, x=tt, deriv=1)
p2 <- predict(lof, x=tt, deriv=2)
plot(p2, type="l"); abline(h=0, lty=3) # not satisfactory:
lof2 <- lokerns(x,y, deriv=2)
plot(lof2, main=
   "lokerns(*, deriv=2) -- much more smooth than predict(*,deriv=2)")
lines(p2); abline(h=0, lty=3)
}
\keyword{utilities}
