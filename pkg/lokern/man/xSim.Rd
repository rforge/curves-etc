\name{xSim}
\alias{xSim}
\non_function{}
\title{Simulated Linear plus Exponential Peak}
\usage{data(xSim)}
\format{
  A vector of 75 numbers between -3.1323 and 4.4505, all rounded to 4
  digits after the decimal.
}
\description{
  This is simulated data, a linear plus an exponential peak.
  In similar form, data like this appears in the smoothing literature
  since at least the eighties.
}
\source{\url{http://www.unizh.ch/biostat/Software/kernf77.html}.
}
\seealso{The example in \code{\link{glkerns}} replicates the
  computations and plots from the source given.
}
\examples{
data(xSim)
plot(xSim, main = "`xSim' - N=75 simulated linear + peak")
}
\keyword{datasets}
