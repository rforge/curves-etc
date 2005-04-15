\name{lokerns}
\alias{lokerns}
\title{Kernel Regression Smoothing with Local Plug-in Bandwidth}
\description{
  Nonparametric estimation of regression functions and their derivatives
  with kernel regression estimators and automatically adapted \bold{local}
  plug-in bandwidth function.
}
\usage{
lokerns(x , y, deriv = 0, n.out=300, x.out=NULL, korder= deriv + 2,
        hetero=FALSE, is.rand=TRUE,
        inputb= is.numeric(bandwidth) && bandwidth > 0,
        m1 = 400, xl=NULL, xu=NULL,
        s=NULL, sig=NULL, bandwidth=NULL)
}
%% FIXME : The argument list of lokerns() and glkerns() are identical
%% -----   and the \arguments{..} description is 99% identical
%% ====> Really should  merge these help pages or make one of them very short!
\arguments{
  \item{x}{vector of design points, not necessarily ordered.}
  \item{y}{vector of observations of the same length as x.}
  \item{deriv}{order of derivative of the regression function to be
    estimated.  Only deriv=0,1,2 are allowed for automatic smoothing,
    whereas deriv=0,1,2,3,4  is possible when smoothing with an input
    bandwidth array. The default value is deriv=0.}
  \item{n.out}{number of output design points where the function has to
    be estimated; default is \code{n.out=300}.}
  \item{x.out}{vector of output design points where the function has to
    be estimated.  The default is an equidistant grid of n.out points
    from min(x) to max(x).}
  \item{korder}{
    nonnegative integer giving the kernel order; it defaults to
    \code{korder = deriv+2} or \eqn{k = \nu + 2} where \eqn{k - \nu}
    must be even.  The maximal possible values are for automatic
    smoothing, \eqn{k \le 4}{k <= 4}, whereas for smoothing with input
    bandwidth array, \eqn{k \le 6}{k <= 6}.
  }
  \item{hetero}{
    logical: if TRUE, heteroscedastic error variables are assumed for
    variance estimation, if FALSE the variance estimation is optimized for
    homoscedasticity. Default value is hetero=FALSE.
  }
  \item{is.rand}{
    logical: if \code{TRUE} (default), random x are assumed and the
    s-array of the convolution estimator is computed as smoothed
    quantile estimators in order to adapt this variability.  If FALSE,
    the s-array is choosen as mid-point sequences as the classical
    Gasser-Mueller estimator, this will be better for equidistant and
    fixed design.
  }
  \item{inputb}{
    logical: if true, a local input bandwidth array is used; if
    \code{FALSE} (default), a data-adaptive local plug-in bandwidths
    array is calculated and used.
  }
  \item{m1}{
    integer, the number of grid points for integral approximation when
    estimating the plug-in bandwidth. The default, 400, may be increased
    if a very large number of observations are available.
  }
  \item{xl, xu}{
    numeric (scalars), the lower and upper bounds for integral
    approximation and variance estimation when estimating the plug-in
    bandwidth. By default (when \code{xl} and \code{xu} are not specified),
    the 87\% middle part of \eqn{[xmin,xmax]} is used.
  }
  \item{s}{
    s-array of the convolution kernel estimator. If it is not given by input
    it is calculated as midpoint-sequence of the ordered design points for
    \code{is.rand=FALSE} or as quantiles estimators of the design density
    for \code{is.rand=TRUE}.
  }
  \item{sig}{variance of the error variables.  If it is not given by
    input or if \code{hetero=TRUE} (no default) it is calculated by a
    nonparametric variance estimator.}
  \item{bandwidth}{
    \emph{local} bandwidth array for kernel regression estimation.  If it is
    not given by input or if \code{inputb=FALSE} a data-adaptive local
    plug-in bandwidth array is used instead.
  }
}
\value{
  a list including used parameters and estimator.
  \item{x}{vector of ordered design points.}
  \item{y}{vector of observations ordered with respect to x.}
  \item{bandwidth}{local bandwidth array which was used for kernel
    regression estimation.}
  \item{x.out}{vector of ordered output design points.}
  \item{est}{vector of estimated regression function or its derivative.}
  \item{sig}{variance estimation which was used for calculating the
    plug-in bandwidths if hetero=TRUE (default) and either inputb=FALSE
    (default) or is.rand=TRUE (default).}
  \item{deriv}{derivative of the regression function which was estimated.}
  \item{korder}{order of the kernel function which was used.}
  \item{xl}{lower bound for integral approximation and variance estimation.}
  \item{xu}{upper bound for integral approximation and variance estimation.}
  \item{s}{vector of midpoint values used for the convolution kernel
    regression estimator.}
}
\details{
  This function calls an efficient and fast algorithm for automatically
  adaptive nonparametric regression estimation with a kernel method.

  Roughly spoken, the method performs a local averaging of the
  observations when estimating the regression function. Analogously, one
  can estimate derivatives of small order of the regression function.
  Crucial for the kernel regression estimation used here is the choice
  the local bandwidth array. Too small bandwidths will lead to a wiggly
  curve, too large ones will smooth away important details.  The
  function lokerns calculates an estimator of the regression function or
  derivatives of the regression function with an automatically chosen
  local plugin bandwidth function. It is also possible to use a local
  bandwidth array which are specified by the user.

  Main ideas of the plugin method are to estimate the optimal bandwidths
  by estimating the asymptotically optimal mean squared error optimal
  bandwidths. Therefore, one has to estimate the variance for
  homoscedastic error variables and a functional of a smooth variance
  function for heteroscedastic error variables, respectively. Also, one
  has to estimate an integral functional of the squared \eqn{k}-th derivative
  of the regression function (\eqn{k=\code{korder}}) for the global
  bandwidth and the squared \eqn{k}-th derivative itself for the local
  bandwidths.

  Here, a further kernel estimator for this derivative is used with a
  bandwidth which is adapted iteratively to the regression function.  A
  convolution form of the kernel estimator for the regression function
  and its derivatives is used. Thereby one can adapt the s-array for
  random design. Using this estimator leads to an asymptotically minimax
  efficient estimator for fixed and random design.  Polynomial kernels
  and boundary kernels are used with a fast and stable updating
  algorithm for kernel regression estimation.

  More details can be found in the references and on
  \url{http://www.unizh.ch/biostat/Software/kernsplus.html}.
}
\references{
  All the references in \code{\link{glkerns}}.
}
\seealso{\code{\link{glkerns}} for \bold{gl}obal bandwidth computation.}
\examples{
data(cars)
lofit <- lokerns(cars$ speed, cars$ dist)
(sb <- summary(lofit$bandwidth))
op <- par(fg = "gray90", tcl = -0.2, mgp = c(3,.5,0))
plot(lofit$band, ylim=c(0,3*sb["Max."]), type="h",#col="gray90",
     ann = FALSE, axes = FALSE)

boxplot(lofit$bandwidth, add = TRUE, at = 304, boxwex = 8,
        col = "gray90",border="gray", pars = list(axes = FALSE))
axis(4, at = c(0,pretty(sb)), col.axis = "gray")
par(op)
par(new=TRUE)
plot(dist ~ speed, data = cars,
     main = "Local Plug-In Bandwidth Vector")
lines(lofit$x.out, lofit$est, col=4)
mtext(paste("bandwidth in [",
            paste(format(sb[c(1,6)], dig = 3),collapse=","),
            "];  Median b.w.=",formatC(sb["Median"])))
}
\keyword{smooth}
