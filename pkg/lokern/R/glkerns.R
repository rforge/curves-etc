### glkerns   kernel regression smoothing with bandwidth selection

glkerns <- function(x, y, deriv = 0, n.out = 300, x.out = NULL,
                    korder = deriv + 2,
                    ihetero = FALSE, irnd = TRUE, inputb = FALSE,
                    m1 = 400, xl = NULL, xu = NULL,
                    s = NULL, sig = NULL, bandwidth = NULL)
{
    ## control and sort x & y (inputgrid and data)
    n <- length(x)
    if (length(y) != n)
        stop("Input grid `x' and data `y' must have the same length.")
    sorvec <- sort.list(x)
    x <- x[sorvec]
    y <- y[sorvec]

    ## compute/sort outputgrid `x.out' (n.out : length of outputgrid)

    if (is.null(x.out)) { 
        n.out <- as.integer(n.out)
        x.out <- seq(min(x), max(x), length = n.out)
    }
    else
        n.out <- length(x.out <- sort(x.out))

    ## korder         kernel order
    ## deriv          derivative of regression function to be estimated
    if (is.null(korder)) korder <- deriv+2
    if (deriv < 0 || deriv > 4)
        stop("Order of derivative `deriv' must be in {0,1,..,4}.")
    if (korder > 6) {
        warning("Kernel order `korder' must be <= 6; set to deriv + 2")
        korder <- deriv+2
    }

    ## ihetero 	homo- or heteroszedasticity of error variables

    ## irnd     	random or non-random t-grid

    ## m1 : discretization for integral functional estimation
    if ((m1 <- as.integer(m1)) < 10)
        stop("number of discretizations `m1' is too small")

    ## xl, xu: lower/upper bound for integral approximation and
    ##		variance estimation
    if (is.null(xl) || is.null(xu)) {
        xl <- 1
        xu <- 0
    }
    
    ## s mid-point grid
    if (is.null(s) || length(s) != n+1)
        s <- as.double(rep(0, n+1))

    ## sig          input variance
    if (is.null(sig)) sig <- as.double(0)

    ## bandwidth    input bandwidth
    if (is.null(bandwidth)) 
    {
        inputb <- as.integer(0)
        bandwidth <- as.double(0)
    } else if (bandwidth <= 0)
        inputb <- as.integer(0)

    if (deriv > 2 && as.integer(inputb) == 0)
        stop("Order of derivative is too large.")
    if (korder > 4 && as.integer(inputb) == 0)
        korder <- deriv+2

    ## internal parameters and arrays for fortran routine
    len1 <- as.integer(length(x)+1)
    work1 <- double(len1*5)
    work2 <- double(m1*3)
    est <- double(n.out)
    irnd1 <- as.integer(1-irnd)

    ## calling fortran routine
    res <- .Fortran("glkerns", 
                    x = as.double(x),
                    y = as.double(y),
                    as.integer(n),
                    x.out = as.double(x.out),
                    as.integer(n.out),
                    deriv = as.integer(deriv),
                    korder = as.integer(korder),
                    ihetero = as.integer(ihetero),
                    irnd = irnd1,
                    as.integer(inputb),
                    m1,
                    xl = as.double(xl),
                    xu = as.double(xu),
                    s = as.double(s),
                    sig = as.double(sig),
                    work1,
                    work2,
                    bandwidth = as.double(bandwidth),         
                    est = as.double(est)
                    )
    return(x = x, y = y, bandwidth = res$bandwidth, x.out = x.out,
           est = res$est, sig = res$sig, deriv = deriv, korder = korder,
           xl = res$xl, xu = res$xu, s = res$s)
}
