## lpepa	local polynomials with Epanechnikov weights
##		for regression functions and derivatives

lpepa <- function(x, y, bandwidth,
                  deriv = 0, n.out = 200, x.out = NULL,
                  order = deriv+1, mnew = 100, var = FALSE)
{
    ## x		inputgrid
    ## y		data

    ## control and sort inputgrid and data
    n <- length(x)
    if (length(y) != n)
        stop("Input grid and data must have the same length.")
    sorvec <- sort.list(x)
    x <- x[sorvec]
    y <- y[sorvec]

    ## bandwidth	bandwidth for estimation

    ## deriv	derivative of regression function to be estimated

    ## n.out	length of outputgrid
    ## x.out	outputgrid

    ## compute and control outputgrid
    if (is.null(x.out)) {
        n.out <- as.integer(n.out)
        x.out <- seq(min(x),max(x),length = n.out)
    }
    else {
        n.out <- length(x.out)
    }

    ## compute vector of bandwiths
    if (length(bandwidth) == 1)
        bandwidth <- as.double(rep(bandwidth,n.out))
    if (length(bandwidth) != n.out)
        stop("Length of bandwith is not equal to length of output grid.")

    ## sort outputgrid and bandwidth
    sorvec <- sort.list(x.out)
    x.out <- x.out[sorvec]
    bandwidth <- bandwidth[sorvec]

    ## order	order of local polynomial approximation
    ## check order, deriv
    if (order < 0) stop("Polynomial order is negative.")
    if (deriv > order)
        stop("Order of derivative is larger than polynomial order.")

    ## mnew		force of restart

    ## var		switch for variance estimation
    var <- as.logical(var)

    ## internal parameters and arrays (see code in ../src/lpepa.f)
    leng <- 10
    nmoms <- as.integer(length(x)/leng+1)
    imoms <- integer(nmoms)
    moms <- double(nmoms*4*(2+order+as.integer(var)))

    ## check internal limitations from fortran routine
    if (2 + order > 12)
        stop("Polynomial order exceeds 10.")

    res <- .Fortran("lpepa",
                    x = as.double(x),
                    y = as.double(y),
                    as.integer(n),
                    bandwidth = as.double(bandwidth),
                    deriv = as.integer(deriv),
                    order = as.integer(order),
                    x.out = as.double(x.out),
                    as.integer(n.out),
                    mnew = as.integer(mnew),
                    as.integer(imoms),
                    as.double(moms),
                    est = double(n.out),
                    as.integer(leng),
                    as.integer(nmoms),
                    var = as.integer(var),
                    est.var = double(n.out),
                    PACKAGE = "lpridge", DUP = FALSE)

    list(x = x, y = y, bandwidth = res$bandwidth,
         deriv = deriv, x.out = x.out, order = order,
         mnew = mnew, var = var, est = res$est, est.var = res$est.var)
}

