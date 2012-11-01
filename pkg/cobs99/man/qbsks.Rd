\name{qbsks}
\alias{qbsks}
\title{Quantile B-Spline with Fixed Knots}
\description{
  Compute B-spline coefficients for regression quantile B-spline with
  stepwise knots selection and quantile B-spline with fixed knots
  \bold{regression spline}, using Ng (1996)'s algorithm.
}
\usage{
qbsks(x,y,w,pw, knots,nknots, degree,Tlambda, constraint,
      n.sub = n1000cut(n), equal,smaller, greater,gradient, coef,maxiter,
      trace, n.equal,n.smaller,n.greater,n.gradient,
      nrq,nl1, neqc, nj0, tau,lam,tmin,kmax,lstart,
      ks,mk.flag, knots.add, ic, print.mesg,
      factor, tol.kn = 1e-6, eps = .Machine$double.eps, print.warn)
}
\arguments{
  \item{x}{numeric vector, sorted increasingly, the abscissa values}
  \item{y}{numeric, same length as \code{x}, the observations.}
  \item{w}{numeric vector of weights, same length as \code{x}, as in
    \code{\link{cobs}}.}
  \item{pw}{penalty weights vector ... ...}%% << FIXME, currently unused
  \item{knots}{numeric vector of knots of which \code{nknots} will be used.}
  \item{nknots}{number of \code{knots} to be used.}% FIXME: Say more in \details
  \item{degree}{integer specifying polynomial degree; must be 1 or 2.}
  \item{Tlambda}{(vector of) smoothing parameter(s) \eqn{\lambda}, see
    \code{\link{drqssbc}}.}
  \item{constraint}{string (or empty) specifying the global constraints;
    see \code{\link{cobs}}.}
  \item{n.sub}{integer, not larger than sample size \code{n}; the
    default has \code{n.sub == n} as long as \code{n} is less than 1000.}
  \item{equal,smaller,greater,gradient}{optional constraint
    specifications; simply passed to \code{\link{drqssbc}}.}
  \item{coef,maxiter,trace}{simply passed to \code{\link{drqssbc}}.}
  \item{n.equal, n.smaller, n.greater, n.gradient}{simply passed to
    \code{\link{drqssbc}}.}
  \item{nrq,nl1,neqc}{integers specifying dimensionalities, directly
    passed to \code{\link{drqssbc}}, see there.}
  \item{nj0}{passed to \code{\link{drqssbc}}, see there.}
  \item{tau}{desired quantile level (in interval \eqn{(0,1)}).}
  \item{lam,tmin,kmax,lstart}{passed to \code{\link{drqssbc}}, see there.}
  \item{ks}{(seems unused)}
  \item{mk.flag}{logical; if true, perform first step knots selection.}
  \item{knots.add}{logical, see \code{\link{cobs}}.}
  \item{ic}{information criterion to use, see \code{\link{cobs}}.}
  \item{print.mesg}{an integer indicating how \code{qbsks2()} should
    print message about its current stages.}
  \item{factor}{passed to \code{\link{drqssbc}}, see there.}
  \item{tol.kn}{\dQuote{tolerance} for shifting the outer knots.}
  \item{eps}{tolerance passed to \code{\link{drqssbc}}.}
  \item{print.warn}{flag indicating if and how much warnings and
    information is to be printed; currently just passed to \code{\link{drqssbc}}.}
}
\details{
  This is an auxiliary function for \code{\link{cobs}(*, lambda = 0)},
  possibly interesting on its own.  This documentation is currently sparse; read
  the source code!
}
\value{
  a list with components
  \item{coef}{..}
  \item{fidel}{..}
  \item{k}{dimensionality of model fit.}
  \item{ifl}{integer ``flag''; the return code.}
  \item{icyc}{integer of length 2, see \code{\link{cobs}}.}
  \item{knots}{the vector of inner knots.}
  \item{nknots}{the number of inner knots.}
  \item{nvar}{the number of ``variables'', i.e. unknowns including
    constraints.}
  \item{lambda}{the penalty factor, chosen or given.}
  \item{pseudo.x}{the pseudo design matrix \eqn{X}, as returned from
    \code{\link{drqssbc}}.}
}
\references{
  Ng, P. (1996)
  An Algorithm for Quantile Smoothing Splines,
  \emph{Computational Statistics \& Data Analysis} \bold{22}, 99--118.
}
\author{Pin Ng; this help page: Martin Maechler.}

\seealso{the main function \code{\link{cobs}}; further
  \code{\link{drqssbc}} which is called from \code{qbsks()}.
}
\keyword{smooth}
\keyword{utilities}
