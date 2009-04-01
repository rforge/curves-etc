### lokerns   kernel regression smoothing with local bandwidth selection

lokerns <- function(x, y=NULL, deriv = 0, n.out = 300, x.out = NULL,
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
    if(is.null(bandwidth)) {
        bandwidth <- double(n.out)
        if(inputb) stop("NULL bandwidth must have inputb = FALSE")
    } else if(length(bandwidth) != n.out)
        stop("'bandwidth' must be of length 'n.out', i.e., ", n.out)

    ## deriv          derivative of regression function to be estimated
    ## korder         kernel order
    if (deriv < 0) stop("Order of derivative is negative.")
    if (deriv > 4 || (deriv > 2 && !inputb))
        stop("Order of derivative is too large.")
    if (is.null(korder) || korder > 6 || (korder > 4 && !inputb))
        korder <- deriv+2

    ## calling fortran routine
    res <- .Fortran("lokerns",			# Fortran arg.names :
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
                    m1,
                    xl = as.double(xl),
                    xu = as.double(xu),
                    s = as.double(s),
                    sig = as.double(sig),
                    work1 = double((n+1)*5),
                    work2 = double(3 * m1),
                    work3 = double(n.out),
                    bandwidth = as.double(bandwidth),
                    est = double(n.out),
                    PACKAGE = "lokern"
                    )[-c(1:2, 16:18)]# all but (x,y) & work*
    if(res$korder != korder)
	warning(paste("'korder' set to ", res$korder,", internally"))

    list(x = x, y = y, bandwidth = res$bandwidth, x.out = x.out,
	 est = res$est, sig = res$sig,
	 deriv = res$deriv, korder = res$korder,
	 xl = res$xl, xu = res$xu, s = res$s)
}
