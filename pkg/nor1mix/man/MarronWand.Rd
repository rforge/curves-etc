\name{MarronWand}
\alias{MarronWand}% to use as "link"
%% Note: These are ``datasets'' , but cannot be handled as strictly
%%       as such since they  *require* the nor1mix package.
\alias{MW.nm1}
\alias{MW.nm2}
\alias{MW.nm2.old}
\alias{MW.nm3}
\alias{MW.nm4}
\alias{MW.nm5}
\alias{MW.nm6}
\alias{MW.nm7}
\alias{MW.nm8}
\alias{MW.nm9}
\alias{MW.nm10}
\alias{MW.nm11}
\alias{MW.nm12}
\alias{MW.nm13}
\alias{MW.nm14}
\alias{MW.nm15}
\alias{MW.nm16}
\title{Marron-Wand Densities as 'norMix' Objects}
\description{
  The fifteen density examples used in Marron and Wand (1992)'s
  simulation study have been used in quite a few subsequent studies,
  can all be written as normal mixtures and are provided here for
  convenience and didactical examples of normal mixtures.
  Number 16 has been added by Jansen et al.
}
\usage{
 MW.nm1	 # Gaussian
 MW.nm2	 # Skewed
 MW.nm2.old # Skewed(old)
 MW.nm3	 # Str Skew
 MW.nm4	 # Kurtotic
 MW.nm5	 # Outlier
 MW.nm6	 # Bimodal
 MW.nm7	 # Separated (bimodal)
 MW.nm8	 # Asymmetric Bimodal
 MW.nm9	 # Trimodal
 MW.nm10 # Claw
 MW.nm11 # Double Claw
 MW.nm12 # Asymmetric Claw
 MW.nm13 # Asymm. Double Claw
 MW.nm14 # Smooth   Comb
 MW.nm15 # Discrete Comb
 MW.nm16 # Distant Bimodal
}
\source{
  They have been translated from Steve Marron's Matlab code,
  % till 2009 or so: \url{http://www.stat.unc.edu/postscript/papers/marron/parameters/nmpar.m},
  now at \url{http://www.unc.edu/depts/stat-or/miscellaneous/marron/parameters/nmpar.m},
  however for number 2, the Matlab code had \code{MW.nm2.old}; and I've
  defined \code{MW.nm2} as from the Annals paper; see also the last
  example below.
}
\author{Martin Maechler}
\references{
  Marron, S. and Wand, M. (1992)
  Exact Mean Integrated Squared Error;
  \emph{Annals of Statistcs} \bold{20}, 712--736.

  For number 16,\cr% in the book Wand & Marron(1995) has
  Janssen, Marron, Verb..., Sarle (1995) .... %% FIXME
}
\examples{
MW.nm10
plot(MW.nm14)

## These are defined as norMix() calls in  ../R/zMarrWand-dens.R
nms <- ls(pat="^MW.nm", "package:nor1mix")
nms <- nms[order(as.numeric(substring(nms,6)))]
for(n in nms) {
   cat("\n",n,":\n"); print(get(n, "package:nor1mix"))
}

## Plot all of them:
op <- par(mfrow=c(4,4), mgp = c(1.2, 0.5, 0), tcl = -0.2,
          mar = .1 + c(2,2,2,1), oma = c(0,0,3,0))
for(n in nms[-17]) plot(get(n, "package:nor1mix"))
mtext("The Marron-Wand Densities", outer= TRUE, font= 2, cex= 1.6)

## and their Q-Q-plots (not really fast):
prob <- ppoints(N <- 100)
for(n in nms[-17])
   qqnorm(qnorMix(prob, get(n, "package:nor1mix")), main = n)
mtext("QQ-plots of Marron-Wand Densities", outer = TRUE,
      font = 2, cex = 1.6)
par(op)

## "object" overview:
cbind(sapply(nms, function(n) { o <- get(n)
      sprintf("\%-18s: K =\%2d; rng = [\%3.1f, \%2.1f]",
              attr(o, "name"), nrow(o),
              min(o[,"mu"] - 3*sqrt(o[,"sig2"])),
              max(o[,"mu"] + 3*sqrt(o[,"sig2"])) )
             }))


## Note that Marron-Wand (1992), p.720  give #2 as
MW.nm2
## the parameters of which at first look quite different from
MW.nm2.old
## which has been the definition in the above "Source" Matlab code.
## It's easy to see that mu_{nm2} = -.3 + 1.2   * mu_{paper},
## and correspondigly,   s2_{nm2} =       1.2^2 * s2_{paper}
## such that they are "identical" apart from scale and location:
op. <- par(mfrow=2:1, mgp= c(1.2,0.5,0), tcl= -0.2, mar=.1+c(2,2,2,1))
plot(MW.nm2)
plot(MW.nm2.old)
par(op.)
}
\keyword{datasets}
\keyword{distribution}
