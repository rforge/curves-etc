\name{norMix}
\alias{norMix}
\alias{[.norMix}
\alias{is.norMix}
\alias{m.norMix}
\alias{mean.norMix}
\alias{var.norMix}
\alias{print.norMix}
\title{Mixtures of Univariate Normal Distributions}
\description{
  Objects of class \code{norMix} represent finite mixtures of
  (univariate) normal (aka Gaussian) distributions.  Methods for
  construction, printing, plotting, and basic computations are provided.
}
\usage{
norMix(mu, sig2 = rep(1,m), sigma = rep(1,m),
       w = NULL, name = NULL, long.name = FALSE)

is.norMix(obj)
m.norMix(obj)
var.norMix(x, \dots)
\method{mean}{norMix}(x, \dots)
\method{print}{norMix}(x, \dots)
\method{[}{norMix}(x, i,j, drop=TRUE)
}
\arguments{
  \item{mu}{numeric vector of length \eqn{K}, say, specifying the means
    \eqn{\mu} of the \eqn{K} normal components.}
  \item{sig2}{\bold{deprecated!} numeric vector of length \eqn{K},
    specifying the variances \eqn{\sigma^2} of the \eqn{K} normal
    components.  Do specify \code{sigma} instead!}
  \item{sigma}{numeric vector of length \eqn{K}, specifying the standard deviations
    \eqn{\sigma} of the \eqn{K} normal components.}
  \item{w}{numeric vector of length \eqn{K}, specifying the mixture proportions
    \eqn{\pi_j}{p[j]} of the normal components, \eqn{j = 1,\dots,K}.
    Defaults to equal proportions}
  \item{name}{optional name tag of the result (used for printing).}
  \item{long.name}{logical indicating if the \code{name} attribute
    should use punctuation and hence be slightly larger than by default.}
  \item{obj,x }{an object of class \code{norMix}.}
  \item{i,j,drop}{for indexing, see the generic \code{\link{[}} extractor function.}
  \item{\dots}{further arguments passed to methods.}
}
\details{
  The (one dimensional) normal mixtures, \R objects of class
  \code{"norMix"}, are constructed by \code{norMix} and tested for by
  \code{is.norMix}.  \code{m.norMix()} returns the number of mixture
  components; the \code{mean()} method for \code{\link{class} "norMix"}
  returns the (theoretical / true) mean \eqn{E[X]} and 
  \code{var.norMix()} % (not a method, call the function explicitly!) 
  the true variance \eqn{E[(X- E[X])^2]} where
  \eqn{X \sim \langle\mathit{norm.mixt}\rangle}{X ~ <norm.mixt>}.

  The subsetting aka \dQuote{extract} method (\code{x[i,j]}; for generic
  \code{[})---when called as \code{x[i,]}---will typically return a
  \code{"norMix"} object unless matrix indexing selects only one row in
  which case \code{x[i, , drop=FALSE]} will return the normal mixture
  (of one component only).

  For further methods (density, random number generation, fitting,
  \dots), see below. %% -> seealso
}
\value{
  \code{norMix} returns objects of class \code{"norMix"} which are
  currently implemented as 3-column matrix with column names \code{mu},
  \code{sigma}, and \code{w}, and further attributes.
  The user should rarely need to access the underlying structure
  directly.
}
\author{Martin Maechler}
\note{
  For \emph{estimation} of the parameters of such a normal mixture,
  we provide a smart parametrization and an efficient implementation of
  the direct MLE or also the EM algorithm, see
  \code{\link{norMixMLE}()} which includes \code{norMixEM()}.
}
\seealso{\code{\link{dnorMix}} for the density,
  \code{\link{pnorMix}} for the cumulative distribution
  and the quantile function (\code{qnorMix}), and
  \code{\link{rnorMix}} for random numbers and
  \code{\link{plot.norMix}}, the plot method.

  \code{\link{MarronWand}} has the Marron-Wand densities as normal mixtures.

  \code{\link{norMixMLE}()} and \code{norMixEM()} provide \emph{fitting}
  of univariate normal mixtures to data.
}
\examples{
ex <- norMix(mu = c(1,2,5))# defaults: sigma = 1, equal proportions ('w')
ex
plot(ex, p.comp = TRUE)# looks like a mixture of only 2; 'p.comp' plots components

## The 2nd Marron-Wand example, see also  ?MW.nm2
ex2 <- norMix(name = "#2 Skewed",
                mu = c(0, .5, 13/12),
	     sigma = c(1, 2/3, 5/9),
		 w = c(.2, .2, .6))

m.norMix  (ex2)
mean      (ex2)
var.norMix(ex2)
(e23 <- ex2[2:3,]) # (with re-normalized weights)
stopifnot(is.norMix(e23),
          all.equal(var.norMix(ex2),     719/1080, tol=1e-14),
          all.equal(var.norMix(ex ),      35/9,    tol=1e-14),
          all.equal(var.norMix(ex[2:3,]), 13/4,    tol=1e-14),
          all.equal(var.norMix(e23), 53^2/(12^3*4),tol=1e-14)
)

plot(ex2, log = "y")# maybe "revealing"
}
\keyword{distribution}
