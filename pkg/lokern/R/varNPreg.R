#### varest.R : Nonparametric Variance Estimator
####
#### S/R interface to the resest() Fortran subroutine

#### Copyright © Martin Maechler (2001).
#### This software is distributed under the terms of the GNU GENERAL
#### PUBLIC LICENSE Version 2, June 1991, see the COPYING file from R,
#### or http://www.gnu.org/copyleft/gpl.html

varest <- function(x,y)
{
    ## Purpose: Nonparametric Leave-1-out Residuals and Variance Estimator
    ##	in the model   y[i] = mu(x(i)) + E[i] ,  E[i] ~ (0, sigma^2), i.i.d
    ## -------------------------------------------------------------------------
    ## Arguments: (x,y)
    ## -------------------------------------------------------------------------
    ##     subroutine resest(t,x,n, res,snr,sigma2)
    ##-----------------------------------------------------------------------
    ## purpose:
    ##
    ##       computes one-leave-out residuals for nonparametric estimation
    ##       of residual variance (local linear approximation followed by
    ##       reweighting)
    ##
    ## parameters:
    ##
    ##     input   t(n)      abscissae (ordered: t(i) <= t(i+1))
    ##     input   x(n)      data
    ##     input   n         length of data ( >2 )
    ##     output  res(n)    residuals at t(1),...,t(n)
    ##     output  snr       explained variance of the true curve
    ##     output  sigma2    estimation of sigma^2 (residual variance)

    ## Author: Martin Maechler, Date:  9 Jul 2001, 14:47
    if(2 >= n <- length(x)) stop("n := length(x)  must be at least 3")
    if(is.unsorted(x)) stop("`x' must be ordered increasingly")
    if(n != length(y)) stop("`x' and `y' must have same length")
    .Fortran("resest",
             as.double(x), as.double(y), n,
             res = double(n),
             snr = double(1),
             sigma2 = double(1))[4:6]
}
