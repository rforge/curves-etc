.BG
.FN lpridge
.TL
local polynomial regression fitting with ridging
.DN
nonparametric estimation of regression functions and their derivatives
via local polynomials and local polynomial ridge regression with 
polynomial weight functions
.CS
lpridge(x, y, bandwidth, deriv=0, n.out=200, x.out=NULL,  
	order=NULL, ridge=NULL, weight="epa", mnew=100,
	var=F)
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
.AG ridge
ridging parameter. The default value performs a slight ridging (see DETAILS below). ridge=0 leads to the local polynomial estimator without ridging.
.AG weight
kernel weight function. The default value is weight="epa" for Epanechnikov weights. Other weights are "bi" for biweights (square of "epa") and "tri" for triweights (cube of "epa"). If weight is a vector, it is interpreted as vector of coefficients of the polynomial weight function. Thus, weight="epa" is equivalent to weight=c(1,0,-1).
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
.AG ridge
ridging parameter used.
.AG weight
vector of coefficients of the kernel weight function.
.AG mnew
force to restart the algorithm after mnew updating steps. 
.AG var
logical flag: whether the variance of the estimator was computed.
.AG est
estimator of the derivative of order deriv of the regression function.
.AG est.var
estimator of the variance of est (proportional to residual variance).
.DT
This function uses a fast and stable algorithm for nonparametric estimation of regression functions and their derivatives via local polynomials and local polynomial ridge regression with polynomial weight functions described in the reference SBE&G below. Several parameters described in SBE&G are fixed either in the fortran routine or in the Splus-function. There, you find comments how to change them.
.br
In S&G a bad finite sample behavior of local polynomials for random design was found, and ridging of the estimator was proposed. In this algorithm, we use a ridging matrix corresponding to the smoothness assumption "The squared difference of the derivative of order deriv of the regression function at the point of estimation and the weighted mean of design points is bounded by the residual variance divided by the ridge parameter."
Thus, without any smoothness assumption, ridge=0 would be appropriate, and for a nearly constant derivative of order deriv, a ridge parameter going to infinity behaves well. For equidistant design, ridging influences the estimator only at the boundary. Asymptotically, the influence of any non-increasing ridge parameter vanishes.
So far, our experience with the choice of a ridging parameter is limited. Therefore we have chosen a default value, which performs a slight modification of the local polynomial estimator (with denotations h=bandwidth, d=deriv, and where n_0=length(x)*mean(bandwidth)/diff(range(x)) is a mean number of observations in a smoothing interval):
.br
ridge = 5*sqrt(n_0)*h^(2*d)/((2*d+3)*(2*d+5))
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
.EX
myfit <- lpridge(x,y,b,ridge=0)		# local polynomials
plot(x,y)
lines(myfit$x.out,myfit$est,col=2)
myridge <- lpridge(x,y,b)		# local pol. ridge
lines(myridge$x.out,myridge$est,col=3)

.KW smoothing  
.WR
