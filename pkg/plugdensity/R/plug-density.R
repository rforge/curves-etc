# system("R COMPILE plugin.f")
system("R SHLIB plugin.c")
print( dyn.load("/u/maechler/src/Gasser-et-al/plugin-density/plugin.so"))


plugin.density <- function(x, nout = 201, xout = NULL)
{
    ## Purpose:  Plug-in density estimate (global bandwidth)
    ## -------------------------------------------------------------------------
    ## Arguments: x: data;
    ## 		nout : how many output values -> used if xout is NULL (default)
    ##          xout : explicit output abscissae
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 Mar 1998, 18:55
    ## calls is.sorted() from package "SfS"
    n <- length(x <- sort(x))
    xout <-
        if(is.null(xout)) {
            dx <- diff(rx <- range(x))
            if(dx < sqrt(.Machine$double.eps)) dx <- mean(abs(rx))/1000
            seq(from=rx[1] - dx/10, to=rx[2] + dx/10, length= nout)
        } else if(!is.sorted(xout))
            sort(xout)
    m <- length(xout)
    r <- .C("plugin",
            x = as.double(x), n=n,
            z = xout, m=m,
            f= double(m),
            h= double(1))
    structure(list(x= r$z, y= r$f, call = sys.call(), bw = r$h, n=n),
              class=c("density", "densityEHpi"))# Eva Herrman plug in
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

data(co2)
plot(dco2 <- density(co2), ylim = c(0, 0.03))
(pdco2 <- plugin.density(co2, xout = dco2$x))
lines(pdco2, col = "red")

plot.density(pdco2)

str(pdco2 <- plugin.density(co2))
xo <- pdco2 $x
str(d.co2 <- density(co2, n = length(xo), from=xo[1],to=max(xo),
                     width= 4 * pdco2$bw))
all.equal(d.co2, pdco2[c("x","y")])
# close: [1] "Component y: Mean relative difference: 0.0009028029"

