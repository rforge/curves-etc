## Idea : Change this to only estimate the bandwidth !
## ----> incorporate into base R

plugin.density <- function(x, nout = 201, xout = NULL)
{
    ## Purpose:  Plug-in density estimate (global bandwidth)
    ## -------------------------------------------------------------------------
    ## Arguments: x: data;
    ## 		nout : how many output values -> used if xout is NULL (default)
    ##          xout : explicit output abscissae
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 Mar 1998, 18:55

    n <- length(x <- sort(x))
        if(is.null(xout)) {
            ## R's density() here extends the range depending on bandwidth !
            dx <- diff(rx <- range(x))
            if(dx < sqrt(.Machine$double.eps)) dx <- mean(abs(rx))/1000
            m <- as.integer(nout)
            xout <- seq(from=rx[1] - dx/10, to=rx[2] + dx/10, length= m)
        } else {
            m <- length(xout)
            if(is.unsorted(xout)) xout <- sort(xout)
        }
    r <- .C("plugin",
            x = as.double(x), n=n,
            z = xout, m=m,
            f= double(m),
            h= double(1),
            PACKAGE = "plugdensity")
    structure(list(x= r$z, y= r$f, call = sys.call(), bw = r$h, n=n),
              class=c("densityEHpi", "density"))# Eva Herrman plug in
}

print.densityEHpi <- function(x, digits = getOption("digits"), ...)
{
    cat("EvaHerrmann plugin density estimate\n call :",
        deparse(x$call),"\n n = ", x$n,
        " ;  estimated (Gaussian) bandwidth h = ",
        format(x$bw, digits = digits),"\n")
    str(x[1:2], digits = digits, ...)
    invisible(x)
}
