\name{lokerns}
\alias{lokerns}
\title{Kernel Regression Smoothing with Local Plug-in Bandwidth} 
\description{
  Nonparametric estimation of regression functions and their derivatives
  with kernel regression estimators and automatically adapted local
  plug-in bandwidth function.
}
\usage{
lokerns(x , y, deriv=0, n.out=300, x.out=NULL, korder=NULL,
        ihetero=FALSE, irnd=TRUE, inputb=FALSE, m1=400, xl=NULL, xu=NULL,
        s=NULL, sig=NULL, bandwidth=NULL)
}
\arguments{
  \item{x}{vector of design points, not necessarily ordered.}
  \item{y}{vector of observations of the same length as x.}
  \item{deriv}{order of derivative of the regression function to be
    estimated. Only deriv=0,1,2 are allowed for automatic smoothing,
    whereas deriv=0,1,2,3,4  is possible
    smoothing with input bandwidth array.
    The default value is deriv=0.
  }
  \item{n.out}{
    number of output design points, where the function has to be estimated.
    The default value is n.out=300.
  }
  \item{x.out}{
    vector of output design points, where the function has to be estimated. 
    The default is an equidistant grid of n.out points from min(x) to max(x).
  }
  \item{korder}{
    kernel order.  The default value is korder=deriv+2. Here, korder-deriv has
    to be an even number.  For automatic smoothing, korder has to be
    equal or less than 4, for smoothing with input bandwidth array,
    korder has to be equal or less than 6. 
  }
  \item{ihetero}{
    logical value: if TRUE, heteroscedastic error variables are assumed for 
    variance estimation, if FALSE the variance estimation is optimized for 
    homoscedasticity. Default value is ihetero=FALSE.
  }
  \item{irnd}{
    logical value: if TRUE, random x are assumed and the s-array of the convolution
    estimator is computed as smoothed quantile estimators in order to adapt
    this variability. If FALSE, the s-array
    is choosen as mid-point sequences as the classical Gasser-Mueller estimator,
    this will be better for equidistant and fixed design. Default value is irnd=TRUE.
  }
  \item{inputb}{
    logical value: if TRUE, a local input bandwidth array is used, if 
    FALSE a data-adaptive local
    plug-in bandwidths array is calculated and used. Default value is inputb=FALSE.
  }
  \item{m1}{
    number of grid points for integral approximation when estimating the plug-in
    bandwidth. Default value is m1=400. This value may be increased if a very
    large number of observations are available.
  }
  \item{xl}{
    lower bound for integral approximation and variance estimation when 
    estimating the plug-in bandwidth. If xl and xu are not given by input 
    the 87\% middle part of [xmin,xmax] is used.
  }
  \item{xu}{
    upper bound for integral approximation and variance estimation when 
    estimating the plug-in bandwidth. If xl and xu are not given by input 
    the 87\% middle part of [xmin,xmax] is used.
  }
  \item{s}{
    s-array of the convolution kernel estimator. If it is not given by input
    it is calculated as midpoint-sequence of the ordered design points for
    irnd=FALSE or as quantiles estimators of the design density for irnd=TRUE.
  }
  \item{sig}{
    variance of the error variables. If it is not given by input or if ihetero=TRUE
    (no default) it is calculated by the nonparametric variance estimator.
  }
  \item{bandwidth}{
    local bandwidth array for kernel regression estimation. If it is not given by
    input or if inputb=FALSE a data-adaptive local plug-in bandwidth array is used 
    instead.
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
    plug-in bandwidths if ihetero=TRUE (default) and either inputb=FALSE
    (default) or irnd=TRUE (default).}
  \item{deriv}{derivative of the regression function which was estimated.}
  \item{korder}{order of the kernel function which was used.}
  \item{xl}{lower bound for integral approximation and variance estimation.}
  \item{xu}{upper bound for integral approximation and variance estimation.}
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
  has to estimate an integral functional of the squared k-th derivative
  of the regression function (k=korder) for the global bandwidth and the
  squared k-th derivative itself for the local bandwidths.

  Here, a further kernel estimator for this derivative is used with a
  bandwidth which is adapted iteratively to the regression function.  A
  convolution form of the kernel estimator for the regression function
  and its derivatives is used. Thereby one can adapt the s-array for
  random design. Using this estimator leads to an asymptotically minimax
  efficient estimator for fixed and random design.  Polynomial kernels
  and boundary kernels are used with a fast and stable updating
  algorithm for kernel regression estimation.

  More details can be found in the refered papers and on the www-page
  http://www.unizh.ch/biostat/Software/kernsplus.html.
}
\references{
  All the references in \code{\link{glkerns}}, and additionally,
  about \emph{local plug-in bandwidth estimators}:
  \cr
  B. Seifert, M. Brockmann, J. Engel, and T. Gasser (1994)
  Fast algorithms for nonparametric curve estimation.
  \emph{J. Computational and Graphical Statistics} \bold{3}, 192--213.
}
\seealso{\code{\link{glkerns}} for \bold{gl}obal bandwidth computation.}
\examples{
data(cars)
attach(cars)
myfit <- lokerns(speed, dist)		
plot(speed, dist)
lines(myfit$x.out, myfit$est, col=2)
}
\keyword{smoothing}
% Converted by Sd2Rd version 1.19.


