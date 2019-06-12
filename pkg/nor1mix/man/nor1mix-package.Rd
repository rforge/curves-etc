\name{nor1mix-package}
\alias{nor1mix-package}
\alias{nor1mix}
\docType{package}
\title{\packageTitle{nor1mix}}% ../DESCRIPTION
\description{
\packageDescription{nor1mix}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{nor1mix}
\packageIndices{nor1mix}

  Note that direct Maximum Likelihood ML (via \code{\link{optim}()}) is typically much
  faster converging (and more reliably detecting convergence correctly),
  notably thanks to a smart re-parametrization: use \code{\link{norMixMLE}()}.
}
\author{
\packageAuthor{nor1mix}

Maintainer: \packageMaintainer{nor1mix}
}
%% \references{
%% ~~ Literature or other references for background information ~~
%% }
\seealso{
  The Marron-Wand examples of normal (gaussian) mixtures \code{\link{MarronWand}}.

  Multivariate distributions from copulas \code{\link[copula]{Mvdc}}
  from the \pkg{copula} package can use \code{\link{norMix}} marginals.
}
\examples{
   example(dnorMix)
   %% TODO
}
\keyword{ package }
