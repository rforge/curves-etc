\name{glkerns}
\alias{glkerns}
\title{Kernel Regression Smoothing with Adaptive Plug-in Bandwidth}
\description{
  Nonparametric estimation of regression functions and their derivatives
  with kernel regression estimators and automatically adapted
  (\bold{gl}obal) plug-in bandwidth.
}
\usage{
glkerns(x, y, deriv = 0, n.out = 300, x.out=NULL, korder=NULL,
        ihetero=FALSE, irnd=TRUE, inputb=FALSE, m1=400, xl=NULL, xu=NULL,
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
    whereas deriv=0,1,2,3,4 is possible when smoothing with a global input
    bandwidth.  The default value is deriv=0.}
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
  \item{ihetero}{logical: if TRUE, heteroscedastic error variables are
    assumed for variance estimation, if FALSE the variance estimation is
    optimized for homoscedasticity.  Default value is ihetero=FALSE.
  }
  \item{irnd}{
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
    \code{irnd=FALSE} or as quantiles estimators of the design density
    for \code{irnd=TRUE}.
  }
  \item{sig}{variance of the error variables.  If it is not given by
    input or if \code{ihetero=TRUE} (no default) it is calculated by a
    nonparametric variance estimator.}
  \item{bandwidth}{
    \emph{global} bandwidth for kernel regression estimation.  If it is
    not given by input or if \code{inputb=FALSE} a data-adaptive global
    plug-in bandwidth is used instead.}
}

\value{
  a list including used parameters and estimator.
  \item{x}{vector of ordered design points.}
  \item{y}{vector of observations ordered with respect to x.}
  \item{bandwidth}{bandwidth which was used for kernel regression estimation.}
  \item{x.out}{vector of ordered output design points.}
  \item{est}{vector of estimated regression function or its derivative.}
  \item{sig}{variance estimation which was used for calculating the
    plug-in bandwidth}
  \item{deriv}{derivative of the regression function which was estimated.}
  \item{korder}{order of the kernel function which was used.}
  \item{xl}{lower bound for integral approximation and variance estimation.}
  \item{xu }{upper bound for integral approximation and variance estimation.  }
  \item{s}{vector of midpoint values used for the convolution kernel
    regression estimator.}
}
\details{
  This function calls an efficient
  and fast algorithm for automatically adaptive nonparametric
  regression estimation with a kernel method.

  Roughly spoken, the method performs a local averaging of the
  observations when estimating the regression function. Analogously, one
  can estimate derivatives of small order of the regression function.
  Crucial for the kernel regression estimation used here is the choice
  of a global bandwidth. Too small bandwidths will lead to a wiggly
  curve, too large ones will smooth away important details.  The
  function glkerns calculates an estimator of the regression function or
  derivatives of the regression function with an automatically chosen
  global plugin bandwidth. It is also possible to use global bandwidths
  which are specified by the user.

  Main ideas of the plugin method are to estimate the optimal bandwidths
  by estimating the asymptotically optimal mean integrated squared error
  optimal bandwidths. Therefore, one has to estimate the variance for
  homoscedastic error variables and a functional of a smooth variance
  function for heteroscedastic error variables, respectively. Also, one
  has to estimate an integral functional of the squared k-th derivative
  of the regression function (k=korder) for the global bandwidth.

  Here, a further kernel estimator for this derivative is used with a
  bandwidth which is adapted iteratively to the regression function.  A
  convolution form of the kernel estimator for the regression function
  and its derivatives is used. Thereby one can adapt the s-array for
  random design. Using this estimator leads to an asymptotically minimax
  efficient estimator for fixed and random design.  Polynomial kernels
  and boundary kernels are used with a fast and stable updating
  algorithm for kernel regression estimation.  More details can be found
  in the references and on
  \url{http://www.unizh.ch/biostat/Software/kernsplus.html}.
}
\references{
  - global plug-in bandwidth estimator:\cr
  T. Gasser, A. Kneip & W. Koehler (1991)
  A flexible and fast method for automatic smoothing.
  \emph{Journal of the American Statistical Association} \bold{86}, 643--652.

  - variance estimation:\cr
  T. Gasser, L. Sroka & C. Jennen-Steinmetz (1986)
  Residual and residual pattern in nonlinear regression.
  \emph{Biometrika} \bold{73}, 625--633.

  - adapting heteroscedasticity:\cr
  E. Herrmann (1997)
  Local bandwidth choice in kernel regression estimation.
  \emph{Journal of Graphical and Computational Statistics} \bold{6}, 35--54.

  - fast algorithm for kernel regression estimator:\cr
  T. Gasser & A. Kneip (1989)
  discussion of Buja, A., Hastie, TRUE. and Tibshirani, R.: Linear smoothers
  and additive models, \emph{The Annals of Statistics} \bold{17}, 532--535.

  B. Seifert, M. Brockmann, J. Engel & T. Gasser (1994)
  Fast algorithms for nonparametric curve estimation.
  \emph{J. Computational and Graphical Statistics} \bold{3}, 192--213.
  
  - on the special kernel estimator for random design point:\cr
  E. Herrmann (1996)
  \emph{On the convolution type kernel regression estimator};
  Preprint 1833, FB Mathematik, Technische Universitaet Darmstadt
  (available from
  \url{http://wwwbib.mathematik.tu-darmstadt.de/Math-Net/Preprints/Listen/shadow/pp1833.html})
}
\seealso{\code{\link{lokerns}} for \bold{lo}cal bandwidth computation.}
\examples{
## Example from kernf77.html
##  (simulated data : linear plus an exponential peak)
x <- c(1.9666, 1.9000, 1.6449, 1.4275, 2.4000, 1.8487,
       1.3383, 1.9514, 1.5722, 1.1978, 1.3109, 1.6940, 1.0816,
       1.1070, 1.2650, 1.1143,  .5717,  .9492,  .3116,  .7942,
        .5670,  .3555,  .3243,  .3242,  .7647,  .4167,  .7583,
        .9462, 1.4337, 1.8820, 2.3054, 2.1040, 3.8500, 3.7052,
       4.3047, 4.4505, 4.3425, 4.3516, 3.7578, 4.0112, 3.3286,
       2.9024, 2.4331, 1.5544, 1.1867,  .6870, -.3630,  .0436,
       -.8197, -.2883,-1.0559, -.9795,-1.7802,-2.2206,-1.3922,
      -1.6511,-1.7694,-1.1309,-2.1912,-1.5785,-2.6189,-1.8125,
      -2.6155,-1.4585,-2.0951,-2.1428,-2.4827,-2.4171,-2.6610,
      -2.6509,-2.6475,-3.1040,-2.7631,-2.9486,-3.1323)
n <- length(x)
tt <- ((1:n) - 1/2)/n # equidistant x
str(gk <- glkerns(tt, x))
## local bandwidth almost identical :
str(lk <- lokerns(tt, x))

cols <- c(gl="PaleGreen", lo="Firebrick")
plot(lk$x.out, lk$bandwidth, axes = F, xlab="", ylab="",
     ylim=c(0,max(lk$bandwidth)), type="h", col = "gray90")
axis(4); mtext("bandwidth(s)", side=4)
lines(lk$x.out, lk$bandwidth, col = cols["lo"], lty = 3)
abline(     h = gk$bandwidth, col = cols["gl"], lty = 4)
par(new=T)
plot(tt, x, main = "global and local bandwidth kernel regression")
lines(gk$x.out, gk$est, col = cols["gl"], lwd = 1.5)
lines(lk$x.out, lk$est, col = cols["lo"])
# the red curve (local bw) is very slightly better
legend(0.7,4.4, c("global bw","local bw"), col = cols, lwd=1)

data(cars)
plot(dist ~ speed, data = cars,
     main = "Global Plug-In Bandwidth")
myfit <- glkerns(cars$ speed, cars$ dist)		
lines(myfit$x.out, myfit$est, col=2)
mtext(paste("bandwidth = ", format(myfit$bandwidth, dig = 4)))
}
\keyword{smooth}
