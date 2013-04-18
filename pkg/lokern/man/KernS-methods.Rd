\name{KernS-methods}
\title{Methods for ("KernS" classed) Results of lokerns() and glkerns()}
\alias{fitted.KernS}
\alias{plot.KernS}
\alias{predict.KernS}
\alias{print.KernS}
\alias{residuals.KernS}
\description{
  Methods for results of \code{\link{glkerns}()} and
  \code{\link{lokerns}()} which are of (S3) class \code{"KernS"}.
}
\usage{
\S3method{fitted}{KernS} (object, \dots)
\S3method{plot}{KernS} (x, \dots)
\S3method{predict}{KernS} (object, x, deriv = 0, trace.lev = 0, \dots)
\S3method{print}{KernS} (x, digits = getOption("digits"), \dots)
\S3method{residuals}{KernS} (object, \dots)
}
\arguments{
  \item{x,object}{an \R object, of S3 class \code{"KernS"}, typically
    result either from \code{\link{glkerns}()} or \code{\link{lokerns}()}.}
  \item{deriv}{integer, \eqn{\ge 0}{>= 0}, specifiying order of
    derivative that should be predicted.}
  \item{digits}{number of significant digits, see \code{\link{print}}.}
  \item{trace.lev}{integer; level of tracing of Fortran level
    computations; see \code{\link{lokerns}}.}
  \item{\dots}{potentially further arguments passed to and from
    methods.}
}
% \details{

% }
\value{
  (differing, depending on the generic function)%% FIXME
}
\seealso{
 \code{\link{glkerns}}, \code{\link{lokerns}}.
}
\examples{
## see ?glkerns and ?lokerns  for now
%% FIXME
}
\keyword{utilities}
