!S COMPILE plugin.f
system("R COMPILE plugin.f")
system("R COMPILE plugin.c")
#-> f77 -O2 -c plugin.f
print( dyn.load("/u/maechler/src/Gasser-et-al/plugin-density/plugin.o"))
# "plugin_"

plug.co2 <- .Fortran("plugin",
                     x=co2, n=length(co2),
                     z=seq(310,360, length=201), m=201:201,
                     f= double(201),
                     h= double(1))
##> neg**non-integral: DOMAIN error

plug.co2 <- .Fortran("plugin",
                     x = sort(co2), n = length(co2),
                     z = seq(310,360, length=201), m = 201:201,
                     f= double(201),
                     h= double(1))
str(plug.co2)
str(data.matrix(plug.co2[c("z","f")]))

plugin.density <- function(x, nout = 201, xout = NULL)
{
  ## Purpose:  Plug-in density estimate (global bandwidth)
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 16 Mar 98, 18:55
  n <- length(x <- sort(x))
  if(is.null(xout)) {
    dx <- diff(rx <- range(x))
    if(dx < .Machine$single.eps) dx <- mean(abs(rx))/1000
    xout <- seq(from=rx[1] - dx/10, to=rx[2] + dx/10, length= nout)
  } else if(!is.sorted(xout))
    xout <- sort(xout)
  m <- length(xout)
  r <- .Fortran("plugin",
                x = as.double(x), n=n,
                z = xout, m=m,
                f= double(m),
                h= double(1))
  structure(list(x= r$z, y= r$f, call = sys.call(), h=r$h, n=n),
            class="densityEHpi")# Eva Herrman plug in
}


plot.densityEHpi <- function(denso, type="l", main= deparse(cal),
                             xlab = deparse(cal[[2]]), ylab = "f(.)", ...)
{
  ## Purpose:  plot method for densityEHpi objects
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 16 Mar 98, 19:05
  cal <- denso$call
  plot.default(denso, type=type, main=main, xlab=xlab, ylab=ylab, ...)
}
plot(plug.co2)

str(pi.co2 <- plugin.density(co2))
xo <- pi.co2 $x
str(d.co2 <- density(co2, n = length(xo), from=xo[1],to=max(xo),
                     width= 4 * pi.co2$h))
all.equal(d.co2, pi.co2[c("x","y")])
# close: [1] "Component y: Mean relative difference: 0.0009028029"

