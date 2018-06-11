.BG
.FN lpepa
.TL
local polynomial regression fitting with Epanechnikov weights
.DN
nonparametric estimation of regression functions and their derivatives
via local polynomials with Epanechnikov weight function
.CS
lpepa(x, y, bandwidth, deriv=0, n.out=200, x.out=NULL,
	order=NULL, mnew=100, var=F)
.RA
.AG x
vector of design points, not necessarily ordered.
.AG y
vector of observations of the same length as x.
.AG bandwidth
bandwidth for nonparametric estimation. Either a number or a vector of the same length as x.out.
.OA
.AG deriv
order of derivative of the regression function to be estimated. The default value is deriv=0.
.AG n.out
number of output design points, where the function has to be estimated.
The default value is n.out=200.
.AG x.out
vector of output design points, where the function has to be estimated. The default value is an equidistant grid of n.out points from min(x) to max(x).
.AG order
order of the polynomial used for local polynomials. The default value is order=deriv+1.
.AG mnew
force to restart the algorithm after mnew updating steps. The default value is mnew=100. For mnew=1 you get a numerically "super-stable" algorithm (see ref. SBE&G below).
.AG var
logical flag: if TRUE, the variance of the estimator proportional to the residual variance is computed (see DETAILS below).
.RT
a list including used parameters and estimator.
.AG x
vector of ordered design points.
.AG y
vector of observations ordered according to x.
.AG bandwidth
vector of bandwidths actually used for nonparametric estimation.
.AG deriv
order of derivative of the regression function estimated. 
.AG x.out
vector of ordered output design points.
.AG order
order of the polynomial used for local polynomials.
.AG mnew
force to restart the algorithm after mnew updating steps. 
.AG var
logical flag: whether the variance of the estimator was computed.
.AG est
estimator of the derivative of order deriv of the regression function.
.AG est.var
estimator of the variance of est (proportional to residual variance).
.DT
This function uses a fast and stable algorithm for nonparametric estimation of regression functions and their derivatives via local polynomials with Epanechnikov weight function described in the reference SBE&G below. In S&G a bad finite sample behaviour of local polynomials for random design was found.
For practical use we propose local polynomial regression fitting with ridging, as implemented in the function "lpridge".  In lpepa, several parameters described in SBE&G are fixed either in the fortran routine or in the Splus-function. There, you find comments how to change them.
.br
For var=T, the variance of the estimator proportional to the residual variance is computed, i.e. the exact finite sample variance of the regression estimator is 
.br
var(est) = est.var*sigma^2.
.SH REFERENCES
- numerical stability and computational speed:
.br
B. Seifert, M. Brockmann, J. Engel & T. Gasser (1994). Fast algorithms for nonparametric curve estimation. J. Computational and Graphical Statistics 3, 192-213.
.br
- statistical properties:
.br
B. Seifert & T. Gasser (1994). Finite sample variance of local polynomials: Analysis and solutions. submitted (also on anonymous ftp: biostat1.unizh.ch and WWW: http://www.unizh.ch/biostat).
.SA
lpridge
.EX
myfit <- lpepa(x,y,b)		# local polynomials
plot(x,y)
lines(myfit$x.out,myfit$est)

.KW smoothing  
.WR
