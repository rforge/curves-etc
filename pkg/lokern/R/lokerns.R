### lokerns   kernel regression smoothing with local bandwidth selection

## auxiliary, factored out of lokerns()
## -- also called via do.call() from predict.KernS()
.lokerns <- function(x,y,x.out,nobs,n.out,deriv,korder,
		     hetero,is.rand,inputb,
		     m1,xl,xu,s,sig,bandwidth, trace.lev)
{
    ## FIXME: These checks should happen in the *caller*, but since we see valgrind problems...
    L1 <- function(x) length(x) == 1L
    stopifnot(L1(nobs), L1(n.out),
              nobs == length(x), nobs == length(y), n.out == length(x.out),
              L1(deriv), L1(korder), L1(hetero), L1(is.rand), L1(inputb),
              L1(m1), L1(xl), L1(xu), L1(sig), L1(trace.lev), trace.lev >= 0)
    r <- .Fortran(lokern_s,			# Fortran argument names (iff differing):
		  x = as.double(x),		# t
		  y = as.double(y),		# x
		  x.out = as.double(x.out),	# tt
		  est	= double(n.out),	# y
		  nobs = as.integer(nobs),	# n
		  n.out= as.integer(n.out),	# m
		  deriv = as.integer(deriv),	# nue
		  korder = as.integer(korder),	# kord
		  hetero = asIntLogic(hetero),	# hetero
		  is.rand= asIntLogic(is.rand),	# isrand
		  inputb = asIntLogic(inputb),
		  iter = as.integer(m1),        # m1 :   number of plug-in iterations on output
		  xl = as.double(xl),		# tl
		  xu = as.double(xu),		# tu
		  s = as.double(s),
		  sig = as.double(sig),
		  work1 = double((nobs+1)*5),	# wn
		  work2 = double(3 * m1),	# w1
		  work3 = double(n.out),	# wm
		  bandwidth=as.double(bandwidth)# ban   # 20
		  , as.integer(trace.lev)       # trace
		  )[-c(1:2, 17:19, 21L)] # keep all, but (x,y), work*,..
    if(r$korder != korder)
	warning(gettextf("'korder' reset from %d to %d, internally",
			 korder, r$korder))
    if(!is.finite(r$iter)) warning("r$iter = ", r$iter," is not finite..")
    else if(r$iter < 0) r$iter <- NA_integer_
    r
}

lokerns <- function(x, ...)
    UseMethod("lokerns")

## modeled along t.test.formula ( ~/R/D/r-devel/R/src/library/stats/R/t.test.R ):
lokerns.formula <- function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(mf$data, parent.frame())))
        mf$data <- as.data.frame(data)
    mf[[1L]] <- quote(stats::model.frame)
    mf$... <- NULL
    mf <- eval(mf, parent.frame())
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")

    r <- lokerns.default(mf[[-response]], mf[[response]], ...)
    r$call <- cl
    r
}


lokerns.default <- function(x, y=NULL, deriv = 0,
                            n.out = 300, x.out = NULL, x.inOut = TRUE,
                            korder = deriv + 2, hetero = FALSE, is.rand = TRUE,
                            inputb = is.numeric(bandwidth) && bandwidth > 0,
                            m1 = 400, xl = NULL, xu = NULL, s = NULL, sig = NULL,
                            bandwidth = NULL, trace.lev = 0, ...)
{
    ## control and sort input (x,y) - new: allowing only y
    xy <- xy.coords(x,y)
    x <- xy$x
    n <- length(x)
    if (n < 3) stop("must have n >= 3 observations")
    if(length(list(...))) warning("extraneous arguments ignored: ", deparse(list(...)))
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
	if(identical(x.inOut, FALSE)) {
	    x.out <- seq(x[1], x[n], length = n.out)
	}
	else { ## construct x.out containing x[] and more
	    if(identical(x.inOut, TRUE))
		seqXmethod <- "aim"	# cheaper than "interpolate"
	    else {			## x.inOut is a character
		if(!is.character(x.inOut))
		    stop("'x.inOut' must be logical or character")
		seqXmethod <- match.arg(x.inOut, eval(formals(seqXtend)$method))
		x.inOut <- TRUE
	    }
	    n.out <- length(x.out <- seqXtend(x, n.out, method = seqXmethod))
	    ind.x <- match(x, x.out)
	}
    }
    else {
	n.out <- length(x.out <- sort(x.out))
	ind.x <- match(x, x.out)## x[] matching x.out[]:
	## FIXME: approximate matching would be better: findInterval() etc
	if(!missing(x.inOut))
	    warning("'x.inOut' is disregarded as 'x.out' has been specified")
	x.inOut <- all(!is.na(ind.x))
	if(x.inOut)
	    seqXmethod <- NA_character_ # don't need seqXtend()
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
	xu <- 0 # so xu < xl  ==> compiled code 'phase 4' computes (tl, tu)
    }

    ## s	mid-point grid :
    if(is.null(s) || length(s) != n+1) s <- double(n+1)

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

    xinL <- if(x.inOut) list(ind.x = ind.x, seqXmethod = seqXmethod)
    structure(c(xy[c("x","y")], # (x,y) possibly unsorted..
		.lokerns(x=x,y=y,x.out=x.out,nobs=n,n.out=n.out,deriv=deriv,
			 korder=korder,hetero=hetero,is.rand=is.rand,
			 inputb=inputb,m1=m1,xl=xl,xu=xu,
			 s=s,sig=sig,bandwidth=bandwidth,trace.lev=trace.lev),
		xinL,
		list(m1 = m1, isOrd = isOrd, ord = if(!isOrd) ord,
		     x.inOut = x.inOut, call = match.call())),
	      class = c("lokerns", "KernS"))
}

#### FIXME:  does only work when 'x.out' was 'x' originally
#### -----   Need better: by default  x.out should contain x as  x.out[ind.x]
fitted.KernS <- function(object, ...) {
    if(object$x.inOut) {
	fit <- object$est[object$ind.x]
	if(object$isOrd) fit else fit[order(object$ord)]
    }
    else stop("'KernS' fit was done with 'x.out' not including data;",
                "\n hence cannot provide fitted values or residuals")
}
residuals.KernS <- function(object, ...) object$y - fitted(object)

print.KernS <- function(x, digits = getOption("digits"), ...)
{
    if(!is.null(cl <- x$call)) {
	cat("Call:\n")
	dput(cl, control=NULL)
    }

    ## This is too cheap, but for now ...  FIXME
    str(unclass(x)[names(x) != "call"],
        digits=digits)

    invisible(x)
}

## predict() relies on this :
stopifnot(identical(names(formals(.lokerns)),
                    names(formals(.glkerns))))

predict.KernS <- function (object, x, deriv = object[["deriv"]],
			   korder = deriv+2, trace.lev = 0, ...)
{
    if(missing(x) && length(dots <- list(...)) == 1 && names(dots) == "newdata") {
	x <- dots$newdata
	if(is.data.frame(x)) {
	    message("using first column of data.frame as 'x'")
	    x <- x[,1L]
	}
	if(!is.numeric(x))
	    stop("Need numeric 'x' for predict(<(lo|gl)kerns>, *)")
	## else use 'x' !
    }
    else if(getRversion() >= "3.3")
	chkDots(...) # useRs using other arguments -- no longer silently!
    if(deriv == object$deriv) {
	if (missing(x) && object$x.inOut) {
	    return(list(x = object[["x"]], y = {
		fit <- object[["est"]][object[["ind.x"]]]
		if(object$isOrd) fit else fit[order(object[["ord"]])]
	    }))
	}
	else if(!missing(x) &&
		!any(is.na(mx <- match(x, object[["x.out"]]))))
	    ## the x's are all "there"
	    return(list(x = x, y = object[["est"]][mx]))
    }
    ## else non-matching 'x' or "different" deriv :

    cl1 <- class(object)[[1]]
    FUN <- switch(cl1,
		  "lokerns" = .lokerns,
		  "glkerns" = .glkerns,
		  stop("invalid class(.)[1]"))
    nf <- names(formals(FUN))
    args <- object[nf[nf != "trace.lev"]]
    args[["trace.lev"]] <- trace.lev
    args[["inputb"]] <- TRUE # we *do* provide the bandwidths
    args[["deriv"]] <- deriv
    args[["korder"]] <- korder
    if(missing(x)) x <- object[["x"]]
    args[["n.out"]] <- length(x)
    if(cl1 == "lokerns") {
	## The "peculiar thing": 'bandwidth' (vector) must
	## be of the length of 'x.out' ! -- using spline (inter/extra)polation
	args[["bandwidth"]] <-
	    spline(x=object[["x.out"]], y=object[["bandwidth"]],
		   method = "natural", xout = x)$y
    }
    args[["x.out"]] <- x
    list(x = x, y = do.call(FUN, args)[["est"]])
}

plot.KernS <- function (x, type = "l", lwd = 2.5, col = 3, ...) {
    if(x$deriv == 0)## sfsmisc:: data and curve; even residuals:
        plotDS(x$x, yd = x$y, ys = list(x = x$x.out, y = x$est), ...)
    else plot(x$x.out, x$est, type=type, lwd=lwd, col=col, ...)
}

lines.KernS <- function (x, ...) lines(x$x.out, x$est, ...)

