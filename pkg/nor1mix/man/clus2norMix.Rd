\name{clus2norMix}
\alias{clus2norMix}
\title{Transform Clustering / Grouping to Normal Mixture}
\description{
  Simple transformation of a clustering or grouping to a normal mixture
  object (class \code{"norMix"}, see, \code{\link{norMix}}.
}
\usage{
clus2norMix(gr, x, name = deparse(sys.call()))
}
\arguments{
  \item{gr}{a grouping/clustering vector with values in
    \eqn{\{1,..,K\}}; possibly a \code{\link{factor}}.}
  \item{x}{numeric vector of (original) data (of the same length as \code{gr}).}
  \item{name}{name for \code{\link{norMix}()} object; constructed from
    the call by default.}
}
\value{
  A call to \code{\link{norMix}()} with \code{(mu, sig2, w)} set to the
  empirical values of the groups (as defined by \code{\link{split}(x,gr)}.
}
\note{
  Via this function, any simple clustering algorithm (such
  \code{\link[cluster]{pam}}) can be used as simple mixture model
  fitting procedure.
}
\author{Martin Maechler, Dec. 2007}
\seealso{\code{\link{norMix}}; further \code{\link[cluster]{pam}()} (or
  \code{\link{clara}()}) from package \pkg{cluster} for sensible clusterings.}
\examples{
x9 <- rnorMix(500, MW.nm9)
library(cluster)
pxc <- pam(x9, k=3)
plot(pxc, which = 2)# silhouette

(nm.p9 <- clus2norMix(pxc$clustering, x9))
plot(nm.p9, p.norm=FALSE)
lines(MW.nm9, col="thistle")
}
\keyword{cluster}
\keyword{models}
