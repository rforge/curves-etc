### glkerns   kernel regression smoothing with bandwidth selection

glkerns <- function(x, y=NULL, deriv = 0, n.out = 300, x.out = NULL,
		    korder = deriv + 2, hetero = FALSE, is.rand = TRUE,
		    inputb = is.numeric(bandwidth) && bandwidth > 0,
		    m1 = 400, xl = NULL, xu = NULL, s = NULL, sig = NULL,
		    bandwidth = NULL)
{
    ## control and sort input (x,y) - new: allowing only y
    xy <- xy.coords(x,y)
    x <- xy$x
    y <- xy$y
    n <- length(x)
    if (n < 3) stop("must have n >= 3 observations")
    if(is.unsorted(x)) {
	sorvec <- sort.list(x)
	x <- x[sorvec]
	y <- y[sorvec]
    }

    ## compute/sort outputgrid 'x.out' (n.out : length of outputgrid)

    if (is.null(x.out)) {
        n.out <- as.integer(n.out)
        x.out <- seq(min(x), max(x), length = n.out)
    }
    else
        n.out <- length(x.out <- sort(x.out))

    if(n.out == 0) stop("Must have 'n.out' >= 1")

    ## hetero	homo- or heteroszedasticity of error variables
    ## is.rand	random or non-random t-grid
    ## inputb	input bandwidth or estimation of plug-in bandwidth

    ## m1 : discretization for integral functional estimation
    if ((m1 <- as.integer(m1)) < 3)# was "10", but fortran has 3
        stop("number of discretizations 'm1' is too small")

    ## xl, xu: lower/upper bound for integral approximation and
    ##		variance estimation
    if (is.null(xl) || is.null(xu)) {
        xl <- 1
        xu <- 0
    }

    ## s	mid-point grid :
    s <- double(if(is.null(s) || length(s) != n+1)  n+1 else s)

    ## sig          input variance
    if (is.null(sig)) sig <- 0. #-> Fortran takes 0 = "compute default"

    inputb <- as.logical(inputb)
    if (is.null(bandwidth) || bandwidth < 0)
        bandwidth <- 0.
    else {
        bandwidth <- as.double(bandwidth[1])
        if (bandwidth == 0 && inputb)
            stop("bandwidth = 0 must have inputb = FALSE")
    }

    ## deriv          derivative of regression function to be estimated
    ## korder         kernel order
    if (deriv < 0 || deriv > 4)
        stop("Order of derivative 'deriv' must be in {0,1,..,4}.")
    if (deriv > 2 && !inputb)
        stop("Order of derivative must be <= 2  if (! inputb).")
    if (is.null(korder))
        korder <- deriv+2
    else if (korder > 6) {
        warning("Kernel order 'korder' must be <= 6; set to deriv + 2")
        korder <- deriv+2
    } else if (korder > 4 && !inputb) {
        warning("Kernel order must be <= 4 if(!inputb); set to deriv+2")
        korder <- deriv+2
    }

    ## calling fortran routine
    res <- .Fortran("glkerns",			# Fortran arg.names :
                    x = as.double(x),		# t
                    y = as.double(y),		# x
                    x.out = as.double(x.out),	# tt
                    as.integer(n),		# n
                    as.integer(n.out),		# m
                    deriv = as.integer(deriv),  # nue
                    korder = as.integer(korder),# kord
                    hetero = as.logical(hetero),# hetero
                    is.rand= as.logical(is.rand),# isrand
                    inputb  = inputb,		# smo
                    iter = m1, # m1; contains the number of plug-in iterations on output
                    xl = as.double(xl),
                    xu = as.double(xu),
                    s = as.double(s),
                    sig = as.double(sig),
                    work1 = double((n+1)*5),
                    work2 = double(m1*3),
                    bandwidth = bandwidth,
                    est = double(n.out),
                    PACKAGE = "lokern")
    if(res$korder != korder)
	warning(gettextf("'korder' reset from %d to %d, internally",
			 korder, res$korder))
    if(res$iter < 0) res$iter <- NA_integer_

    list(x = x, y = y, bandwidth = res$bandwidth, x.out = x.out,
	 est = res$est, sig = res$sig,
	 deriv = res$deriv, korder = res$korder, iter = res$iter,
	 xl = res$xl, xu = res$xu, s = res$s)
}
