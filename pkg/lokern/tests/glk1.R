require(lokern)
data(xSim)
xSim # `test' for the dataset

n <- length(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x

str(gk <- glkerns(tt, xSim))
summary(gk$est)
gk$bandwidth
glkerns(tt,xSim, deriv = 1)$bandwidth
glkerns(tt,xSim, deriv = 2)$bandwidth

## Q: Can we have *same* kernel, *same* bandwidth  with different  'deriv'
##    similarly to smooth.spline() ?
##
## Answer: not really,  mainly because have not enough choices
##     (nu,k_{ord}), i.e., because currentl,  nu - k  must be even ...

p.3glks <- function(x.dat, y.dat, kord, is.rand=FALSE, useBandwidth, bw.factor = 1.8)
{
    ## Purpose: Plot  glkerns(*,  deriv = {0, 1, 2}
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  2 Jul 2009, 09:24

    if(!missing(useBandwidth) && is.numeric(useBandwidth) && useBandwidth > 0)
        bw <- useBandwidth
    else {
        ## Determine the fixed bandwidth :
        bw0 <- glkerns(x.dat, y.dat, kord=kord, is.rand=is.rand)$bandwidth
        bw <- bw0 * bw.factor     # more smoothing for the derivatives
    }

    ## Estimates for   g, g' , g'' :
    gk0 <- glkerns(x.dat, y.dat, kord=kord, is.rand=is.rand, bandw = bw)
    gk1 <- glkerns(x.dat, y.dat, kord=kord, is.rand=is.rand, bandw = bw, deriv=1)
    ##-> warning kord := 3  ===> indeed: do NOT have internally consistent deriv's
    gk2 <- glkerns(x.dat, y.dat, kord=kord, is.rand=is.rand, bandw = bw, deriv=2)

    ## Plots ---------------

    op <- par(mfrow= c(3,1), mgp = c(1.25, 0.6, 0),
              mar = c(3,3,2.5,1) + .1, oma = c(0,0, 2, 0))
    on.exit(par(op))
    with(gk0, { plot(est ~ x.out, type = "l", main = expression(hat(g)(.)),
                     col=2, lwd = 1.5)
                points(y ~ x, cex = 0.5)
                mtext(substitute(list(bw == B,k[ord] == K),
                                 list(B = formatC(bandwidth), K = korder)),
                      adj = 1, line = .25) })
    with(gk1, { plot(est ~ x.out, type = "l", main = expression(widehat(g*minute)(.)),
                     col=2, lwd = 1.5)
                abline(h = 0, col = "gray", lty=3)
                mtext(substitute(list(bw == B, k[ord] == K),
                                 list(B = formatC(bandwidth), K = korder)),
                      adj = 1, line = .25) })
    with(gk2, { plot(est ~ x.out, type = "l", main = expression(widehat(g*second)(.)),
                     col=2, lwd = 1.5)
                abline(h = 0, col = "gray", lty=3)
                mtext(substitute(list(bw == B, k[ord] == K),
                                 list(B = formatC(bandwidth), K = korder)),
                      adj = 1, line = .25) })
    mtext(sprintf("glkerns(*, deriv = {0, 1, 2}, bandwidth = <fixed>, kord = %d)",
                  kord),
          line = 0.5, outer = TRUE, cex = par("cex.main"), font = par("font.main"))

    invisible(list(gk0 = gk0, gk1 = gk1, gk2 = gk2))
}

p.3glks(tt, xSim, kord = 4)
p.3glks(tt, xSim, kord = 4, useB = 0.15)

p.3glks(tt, xSim, kord = 5, useB = 0.1) ## fixme: ^g uses k =2

p.3glks(tt, xSim, kord = 6, useB = 0.2)

## "FIXME" visually compare with numerical derivatives (e.g. from splines).
