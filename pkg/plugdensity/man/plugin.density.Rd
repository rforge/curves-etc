\name{plugin.density}
\alias{plugin.density}
\alias{print.densityEHpi}
\title{Kernel Density Estimation by Plug-In Bandwidth Selection}
\description{
  The function provides kernel density estimation with iterative plug-in
  bandwidth selection.
}
\usage{
plugin.density(x, nout = 201, xout = NULL)
}
\arguments{
  \item{x}{numeric vector whose density is to be estimated.}
  \item{nout}{ ~~Describe \code{nout} here~~ }
  \item{xout}{ ~~Describe \code{xout} here~~ }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  an object of class \code{"densityEHpi"} inheriting also from class
  \code{"\link{density}"}.
  It is a list with components
  \item{comp1 }{.......}
  \item{... }{.......}
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

plot.density(pdco2)

str(pdco2 <- plugin.density(co2))
xo <- pdco2 $x
str(d.co2 <- density(co2, n = length(xo), from=xo[1],to=max(xo),
                     width= 4 * pdco2$bw))
all.equal(d.co2, pdco2[c("x","y")])
# close: [1] "Component y: Mean relative difference: 0.0009028029"

}
\keyword{distribution}
\keyword{smooth}
