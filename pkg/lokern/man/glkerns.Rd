\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\name{glkerns}
\title{Kernel Regression Smoothing with Adaptive Plug-in Bandwidth}
\alias{glkerns}
\alias{glkerns.default}
\alias{glkerns.formula}
\description{
  Nonparametric estimation of regression functions and their derivatives
  with kernel regression estimators and automatically adapted
  (\bold{gl}obal) plug-in bandwidth.
}
\usage{
glkerns(x, \dots)

\method{glkerns}{default}(x, y=NULL, deriv = 0, n.out = 300, x.out=NULL, x.inOut = TRUE,
        korder= deriv + 2, hetero=FALSE, is.rand=TRUE,
        inputb= is.numeric(bandwidth) && bandwidth > 0,
        m1 = 400, xl=NULL, xu=NULL,
        s=NULL, sig=NULL, bandwidth=NULL, trace.lev = 0, \dots)

\method{glkerns}{formula}(formula, data, subset, na.action, \dots)
}
%% FIXME : The argument list of lokerns() and glkerns() are identical
%% -----   and the \arguments{..} description is 99% identical
%% ====> Really should  merge these help pages or make one of them very short!
\arguments{
  \item{x}{vector of design points, not necessarily ordered.}
  \item{y}{vector of observations of the same length as \code{x}.}
  \item{deriv}{order of derivative of the regression function to be
    estimated.  Only deriv=0,1,2 are allowed for automatic smoothing,
    whereas deriv=0,1,2,3,4 is possible when smoothing with a global input
    bandwidth.  The default value is deriv=0.}
  \item{n.out}{number of output design points where the function has to
    be estimated; default is \code{n.out=300}.}
  \item{x.out}{vector of output design points where the function has to
    be estimated.  The default is an equidistant grid of n.out points
    from min(x) to max(x).}
  \item{x.inOut}{logical or character string indicating if \code{x.out}
    should contain the input \code{x} values.
    Note that this argument did not exist, equivalently to being
    \code{FALSE}, up to \pkg{lokern} version \code{1.0-9}.

    In order for \code{\link{residuals}()} or \code{\link{fitted}()}
    methods to be applicable, it must be \code{TRUE} or a character
    string specifying one of the \code{methods}s of
    \code{\link{seqXtend}} (package \CRANpkg{sfsmisc}). The default,
    \code{TRUE} corresponds to method \code{"aim"}.
  }
  \item{korder}{
    nonnegative integer giving the kernel order \eqn{k}; it defaults to
    \code{korder = deriv+2} or \eqn{k = \nu + 2} where \eqn{k - \nu}
    must be even.  The maximal possible values are for automatic
    smoothing, \eqn{k \le 4}{k <= 4}, whereas for smoothing with input
    bandwidth, \eqn{k \le 6}{k <= 6}.
  }
  \item{hetero}{logical: if TRUE, heteroscedastic error variables are
    assumed for variance estimation, if FALSE the variance estimation is
    optimized for homoscedasticity.  Default value is hetero=FALSE.
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
    \code{FALSE} (by default when \code{bandwidth} is not specified), a
    data-adaptive local plug-in bandwidths array is calculated and used.
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
    input or if \code{hetero=TRUE} it is calculated by a
    nonparametric variance estimator.}
  \item{bandwidth}{
    \emph{global} bandwidth for kernel regression estimation.  If it is
    not given by input or if \code{inputb=FALSE} a data-adaptive global
    plug-in bandwidth is used instead.}
  \item{trace.lev}{integer indicating how much the internal (Fortran
    level) computations should be \dQuote{traced}, i.e., be reported.
    The default, \code{0}, does not print anything.}

%%% formula method:
  \item{formula}{a \code{\link{formula}} of the form \code{y ~ pred},
    specifying the response variable \code{y} and predictor variable
    \code{pred} which must be in \code{data}.}
  \item{data}{an optional matrix or data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{\dots}{for the \code{formula} method: Optional arguments all
    passed to \code{glkerns.default()}.}
}
\value{
  an object of class(es) \code{c("glkerns", "KernS")}, which is
  a list including used parameters and estimator, containing among others
  \item{x}{vector of ordered design points.}
  \item{y}{vector of observations ordered with respect to x.}
  \item{bandwidth}{bandwidth which was used for kernel regression estimation.}
  \item{x.out}{vector of ordered output design points.}
  \item{est}{vector of estimated regression function or its derivative
    (at \code{x.out}).}
  \item{sig}{variance estimation which was used for calculating the
    plug-in bandwidth}
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
  has to estimate an integral functional of the squared \eqn{k}-th derivative
  of the regression function (\eqn{k=\code{korder}}) for the global bandwidth.

  Here, a further kernel estimator for this derivative is used with a
  bandwidth which is adapted iteratively to the regression function.  A
  convolution form of the kernel estimator for the regression function
  and its derivatives is used. Thereby one can adapt the s-array for
  random design. Using this estimator leads to an asymptotically minimax
  efficient estimator for fixed and random design.  Polynomial kernels
  and boundary kernels are used with a fast and stable updating
  algorithm for kernel regression estimation.  More details can be found
  in the references and on
  \url{http://www.biostat.uzh.ch/research/software/kernel.html}
  %was http://www.unizh.ch/biostat/Software/kernsplus.html}.
}
\author{
  The formula method was added after proposals by Andri Signorell.
}
\references{
  - global plug-in bandwidth estimator:\cr
  Theo Gasser, Alois Kneip & Walter Koehler (1991)
  A flexible and fast method for automatic smoothing.
  \emph{Journal of the American Statistical Association} \bold{86}, 643--652.
%  \url{http://www.jstor.org/stable/2290393}

  Muller, H.-G. (1984)
  Smooth optimum kernel estimators of densities, regression curves and modes.
  \emph{The Annals of Statistics} \bold{12}, 766--774.
%  \url{http://www.jstor.org/stable/2241411}

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
  Preprint 1833, FB Mathematik, Technische Universitaet Darmstadt;
  currently available from
  \url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.32.6383}
  %% no longer available from
  %% \url{http://wwwbib.mathematik.tu-darmstadt.de/Math-Net/Preprints/Listen/shadow/pp1833.html})
}
\seealso{\code{\link{lokerns}} for \bold{lo}cal bandwidth computation.
  \code{\link{plot.KernS}} documents all the methods for \code{"KernS"}
  classed objects.

  The \code{\link{demo}} for computing derivatives, \code{demo("glk-derivs")}.
}
\examples{
data(xSim)## linear plus an exponential peak, see help(xSim)
n <- length(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x ==> is.rand = FALSE
gk <- glkerns(tt, xSim, is.rand = FALSE)
gk # print method
plot(gk) # nice plot() method
if(require("sfsmisc")) {
  TA.plot(gk)
} else { plot(residuals(gk) ~ fitted(gk)); abline(h = 0, lty=2) }
qqnorm(residuals(gk), ylab = "residuals(gk)")

cat("glkerns() bandwidth:",format(gk$bandwidth, dig=10),"\n")
## local bandwidth: fit is very similar :
(lk <- lokerns(tt, xSim, is.rand = FALSE))
nobs(lk)

cols <- c(gl="PaleGreen", lo="Firebrick")
plot(lk$x.out, lk$bandwidth, axes = FALSE, xlab="", ylab="",
     ylim=c(0,max(lk$bandwidth)), type="h", col = "gray90")
axis(4); mtext("bandwidth(s)", side=4)
lines(lk$x.out, lk$bandwidth, col = cols["lo"], lty = 3)
abline(     h = gk$bandwidth, col = cols["gl"], lty = 4)
par(new=TRUE)
plot(tt, xSim, main = "global and local bandwidth kernel regression")
lines(gk$x.out, gk$est, col = cols["gl"], lwd = 1.5)
lines(lk$x.out, lk$est, col = cols["lo"])
# the red curve (local bw) is very slightly better
legend(0.7,4.4, c("global bw","local bw"), col = cols, lwd=1)

## This should look
op <- par(mfrow = c(3,1), mar = .1 + c(4,4,2,1), oma = c(0,0,3,0),
          mgp = c(1.5, 0.6,0))
plot(gk, main = expression(paste("Data & ", hat(f))))
## calling extra plot() method
gk1 <- glkerns(tt, xSim, deriv = 1, is.rand = FALSE)
plot(gk1$x.out, gk1$est, col = "green", lwd = 1.5, type = "l",
     main = expression(widehat(paste(f,"'"))))
abline(h=0, col="gray", lty = 3)
gk2 <- glkerns(tt, xSim, deriv = 2, is.rand = FALSE)
plot(gk2$x.out, gk2$est, col = "orange", lwd = 1.5, type = "l",
     main = expression(widehat(paste(f,"''"))))
abline(h=0, col="gray", lty = 3)
mtext("Example from www.unizh.ch/biostat/..../kernf77.html",side=3,
      outer = TRUE, cex = 1, font = par("font.main"))

par(op)
data(cars)
plot(dist ~ speed, data = cars,
     main = "Global Plug-In Bandwidth")
## these two are equivalent
m1glk <- glkerns(dist ~ speed, data = cars)
m.glk <- glkerns(cars$ speed, cars$ dist)
lines(m.glk, col=2) # using the lines() method
mtext(paste("bandwidth = ", format(m.glk$bandwidth, dig = 4)))
ii <- names(m1glk) != "call"
stopifnot(all.equal(m1glk[ii], m.glk[ii], tol = 1e-15))
}
\keyword{smooth}
