\name{plugin.density}
\alias{plugin.density}
\alias{print.densityEHpi}
\title{Kernel Density Estimation by Plug-In Bandwidth Selection}
\description{
  The function provides kernel density estimation with iterative plug-in
  bandwidth selection.
}
\usage{
plugin.density(x, nout = 201, xout = NULL, na.rm = FALSE)
}
\arguments{
  \item{x}{vector of numbser whose density is to be estimated.}
  \item{nout}{integer specifying the number of equispaced \code{xout}
    values to use \emph{only when} \code{xout = NULL} (as by default).}
  \item{xout}{numeric vector of abscissa values at which the density is
    to be evaluated. By default, an equispaced sequence of values
    covering (slightly more than) the range of \code{x}.}
  \item{na.rm}{logical; if \code{TRUE}, missing values are removed
    from \code{x}. If \code{FALSE} any missing values cause an error.}
}
%\details{
% ...
%}
\value{
  an object of class \code{"densityEHpi"} inheriting also from class
  \code{"\link{density}"}.
  It is a list with components
  \item{x}{the \code{n} coordinates of the points where the density is
    estimated.}
  \item{y}{the estimated density values.}
  \item{bw}{the bandwidth used.}
  \item{n}{the sample size after elimination of missing values.}
  \item{call}{the call which produced the result.}
  \item{data.name}{the deparsed name of the \code{x} argument.}
%% 'density()' also gives
%  \item{has.na}{logical, for compatibility (always \code{FALSE}).}
}
\references{
  J. Engel, Eva Herrmann and Theo Gasser (1994).
  An iterative bandwidth selector for kernel estimation of densities and
  their derivatives.
  \emph{Journal of Nonparametric Statistics} \bold{4}, 21--34.
}
\author{Algorithm and C code: Eva Herrmann
  \email{eherrmann@mathematik.tu-darmstadt.de};  R interface:
  Martin Maechler \email{maechler@R-project.org}.
}
\seealso{ \code{\link{density}}. }
\examples{
data(co2)
plot(dco2 <- density(co2), ylim = c(0, 0.03))
(pdco2 <- plugin.density(co2, xout = dco2$x))
lines(pdco2, col = "red")

plot(pdco2)# calls 'plot.density' method

str(pdco2 <- plugin.density(co2))
xo <- pdco2 $x
str(d.co2 <- density(co2, n = length(xo), from=xo[1],to=max(xo),
                     width= 4 * pdco2$bw))
nms <- c("x","y", "bw", "n")
all.equal(d.co2[nms], pdco2[nms])
## are close: "Component 2 (= 'y'): Mean relative difference: 0.0009..."
}
\keyword{distribution}
\keyword{smooth}
