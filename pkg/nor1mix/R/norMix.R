####---- Normal Mixtures  "NorMix" -------
####---- ~~~~~~~~~~~~~~~   ######  -------
#### Object-oriented  S/R - functions for dealing with normal mixtures.
#### -------------------------------------------------------------------------
#### Author: Martin Maechler, 20 Mar 1997
#### -------------------------------------------------------------------------

NorMix <- function(mu, sig2, w = NULL, name = NULL, long.name = FALSE)
{
    ## Purpose: Constructor for 'NorMix' (normal mixture) objects
    ## -------------------------------------------------------------------------
    ## Arguments: mu: vector of means;  sig2: vector of  variances  sigma^2
    ##		w : vector of weights (adding to 1) -- default: equal
    ##		name : name attribute; constructed from (mu,sig2) by default
    ##		long.name : logical used for default \code{name} construction
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 97, 14:58

    if(!is.numeric(mu)) stop("'mu' must be numeric!")
    m <- length(mu)
    if(length(sig2) == 1) sig2 <- rep(sig2, m)
    if(length(sig2) != m || !is.numeric(sig2)|| any(sig2 <=0))
        stop("'sig2' = sigma^2  must be > 0 with same length as 'mu'")
    if(is.null(w))
        w <- rep(1/m, m)
    else {
        if(length(w) != m || !is.numeric(w) || any(w<0))
            stop("'w' must be >= 0  with same length as 'mu'")
        s <- sum(w)
        if(abs(s-1) > 10*.Machine$double.eps) w <- w / sum(w)
    }
    if(is.null(name)) {
        sformat <- function(v) sapply(v, format, digits=1)
        pPar <- function(pp) {
            if(long.name)
                paste("(",paste(sformat(pp),  collapse= ","),")", sep="")
            else
                paste(sformat(pp), collapse= "")
        }
        name <- paste("NM",format(m),".",
                      pPar(mu), "_", pPar(sig2), sep="")
    }
    structure(name = name, class = "NorMix",
              .Data = cbind(mu = mu, sig2 = sig2, w = w))
}

is.NorMix <- function(obj)
{
  ## Purpose: is 'obj' a "NorMix", i.e.  Normal Mixture object ?
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:38
  inherits(obj, "NorMix") &&
  (!is.null(w <- obj[,"w"])) &&
  is.numeric(w) && all(w >= 0) && abs(sum(w)-1) < 1000*.Machine$double.eps
}

m.NorMix <- function(obj) nrow(obj) ##  Number of components of  normal mixture

mean.NorMix <- function(obj)
{
  ## Purpose: Return "true mean", i.e., the expectation of  a normal mixture.
  if(!is.NorMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  obj[,"w"] %*% obj[,"mu"]
}

var.NorMix <- function(obj)
{
  ## Purpose: 'true' Variance, i.e. E[(X- E[X])^2]  for X ~ normal mixture.
  if(!is.NorMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  w <- obj[,"w"]
  mu <- w %*% obj[,"mu"]
  w %*% (obj[,"sig2"] + (obj[,"mu"]-mu)^2)
}

print.NorMix <- function(obj, ...)
{
  ## Purpose: print method for  "NorMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:02
  has.nam <- !is.null(nam <- attr(obj,"name"))
  cat("'Normal Mixture' object",
      if(has.nam) paste("\t ``", nam, "''", sep=''), "\n")
  if(has.nam) attr(obj, "name") <- NULL
  cl <- class(obj);  cl <- cl[ cl != "NorMix"] #- the remaining classes
  class(obj) <- if(length(cl)>0) cl else NULL
  NextMethod("print", ...)
  invisible(obj)
}

dNorMix <- function(obj, x = NULL, xlim = NULL, n = 511,...)
{
  ## Purpose: density evaluation for "NorMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Arguments: obj: Normal Mixture object;  x: abscissa values where to eval
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.NorMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  if(is.null(x)) {
    if(is.null(xlim)) { ##-- construct "reasonable" abscissa values
      xlim <- mean.NorMix(obj) + c(-3,3)*sqrt(var.NorMix(obj))
    }
    x <- seq(xlim[1], xlim[2], length = n)
  }
  m <- m.NorMix(obj) #-- number of components
  y <- numeric(length(x))
  w <- obj[,"w"]; mu <- obj[,"mu"]; sig2 <- obj[,"sig2"]
  for(j in 1:m)
    y <- y + w[j] * dnorm(x, mean = mu[j], sd = sqrt(sig2[j]))
  list(x = x, y = y)
}


plot.NorMix <- function(obj, main = attr(obj,"name"), n = 511, xx = NULL,
			p.norm = TRUE, p.h0 = TRUE,
                        norm.col = 2, h0.col=3, ...)
{
  ## Purpose: plot method for  "NorMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  d.o <- dNorMix(obj, n=n, x = xx); x <- d.o$x
  if(p.norm)
    dn  <- dnorm(x, mean = mean.NorMix(obj), sd = sqrt(var.NorMix(obj)))
  plot(d.o, type='l', ylim = c(0,max(d.o$y, if(p.norm) dn)),
       xlab = 'x', ylab = 'f(x)', main = main, lwd = 1.4, ...)
  if(p.norm) lines(x, dn, lwd=0.1, lty=2, col=norm.col)
  if(p.h0)   abline(h=0,  lwd=0.1, lty=3, col=h0.col)
  invisible(obj)
}

r.NorMix <- function(obj, x = NULL, xlim = NULL, n = 511,
		     xy.return = TRUE, ...)
{
  ## Purpose: Compute r := f / f0; f = normal mixture; f0 = "best" normal approx
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.NorMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  m <- m.NorMix(obj) #-- number of components
  d.o <- dNorMix(obj, x=x, xlim=xlim, n=n)
  dn  <- dnorm(d.o$x, mean = mean.NorMix(obj), sd = sqrt(var.NorMix(obj)))
  if(xy.return) list(x = d.o$x, y= d.o$y / dn, f0 = dn) else d.o$y / dn
}

### ---> NorMix-ex.R  for calling these
