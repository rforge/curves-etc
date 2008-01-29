### lokerns   kernel regression smoothing with local bandwidth selection

if(FALSE) ## FIXME: use  sfsmisc::seqXtend()  -- !
seqX <- function(x, n.out, xlim = range(x))
{
    ## sequence() of desired length 'n.out'  *containing*  x
}

lokerns <- function(x, y=NULL, deriv = 0,
                    n.out = 300, x.out = NULL, x.inOut = TRUE,
		    korder = deriv + 2, hetero = FALSE, is.rand = TRUE,
		    inputb = is.numeric(bandwidth) && bandwidth > 0,
		    m1 = 400, xl = NULL, xu = NULL, s = NULL, sig = NULL,
		    bandwidth = NULL)
{
    ## control and sort input (x,y) - new: allowing only y
    xy <- xy.coords(x,y)
    x <- xy$x
    n <- length(x)
    if (n < 3) stop("must have n >= 3 observations")
    x.isInd <- !is.null(xy$xlab) && xy$xlab == "Index"
    isOrd <- x.isInd || !is.unsorted(x)
    if(isOrd)
        y <- xy$y
    else {
        ord <- sort.list(x)
        x <- x[ord]
	y <- y[ord]
    }

    ## compute/sort outputgrid 'x.out' (n.out : length of outputgrid)

    if (is.null(x.out)) {
        n.out <- as.integer(n.out)
        if(!x.inOut){
            x.out <- seq(x[1], x[n], length = n.out)
        }
        else {# construct x.out containing x[] and more
            nDup <- !duplicated(x)
            nu <- length(ux <- x[nDup])
            ind.x <- cumsum(nDup) # ! ==> x === ux[ind.x]
            if(n.out <= nu) { # take just unique x[]
                x.out <- ux
                n.out <- nu
            } else { # extend ux[] to more points -> x.out[]
                ## too cheap: should adapt to diff(ux) and depend on (nu,n.out)!
                nn <- n.out - nu + 2
                x.out <- c(ux, seq(x[1], x[n], length = nn)[-c(1,nn)])
                ii <- sort.list(x.out)
                x.out <- x.out[ii]
                ind.x <- which(ii <= nu)[ind.x]# ==> x === x.out[ind.x]
            }
        }
    }
    else {
        n.out <- length(x.out <- sort(x.out))
        ind.x <- match(x, x.out)## x[] matching x.out[]:
        ## FIXME: approximate matching would be better: findInterval() etc
        x.inOut <- all(!is.na(ind.x))
    }

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
    s <- if (is.null(s) || length(s) != n+1) double(n+1) else as.double(s)

    ## sig      input variance
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
                    n,				# n
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

    structure(c(xy[c("x","y")], res, # (x,y) possibly unsorted..
                list(isOrd = isOrd, ord = if(!isOrd) ord,
                     x.inOut = x.inOut, ind.x = ind.x,
                     call = match.call())),
              class = c("lokerns", "KernS"))
}

#### FIXME:  does only work when 'x.out' was 'x' originally
#### -----   Need better: by default  x.out should contain x as  x.out[ind.x]
fitted.KernS <- function(object, ...) {
    if(object$x.inOut)
        with(object,
         {
             fit <- est[ind.x]
             if(isOrd) fit else fit[order(ord)]
         })
    else stop("'KernS' fit was done with 'x.out' not including data;",
                "\n hence cannot provide fitted values or residuals")
}
residuals.KernS <- function(object, ...) object$y - fitted(object)

