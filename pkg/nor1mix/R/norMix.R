####-*- mode: R; kept-old-versions: 12;  kept-new-versions: 30; -*-

####---- Normal Mixtures  "norMix" -------
####---- ~~~~~~~~~~~~~~~   ######  -------
#### Object-oriented  S/R - functions for dealing with 1D normal mixtures.
#### -------------------------------------------------------------------------
#### Author: Martin Maechler, 20 Mar 1997
#### -------------------------------------------------------------------------

norMix <- function(mu, sig2 = rep(1, m), w = NULL,
		   name = NULL, long.name = FALSE)
{
    ## Purpose: Constructor for 'norMix' (normal mixture) objects
    ## -------------------------------------------------------------------------
    ## Arguments: mu: vector of means;	sig2: vector of	 variances  sigma^2
    ##		w : vector of weights (adding to 1) -- default: equal
    ##		name : name attribute; constructed from (mu,sig2) by default
    ##		long.name : logical used for default \code{name} construction
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 97, 14:58

    if(!is.numeric(mu)) stop("'mu' must be numeric!")
    m <- length(mu)
    if(length(sig2) == 1) sig2 <- rep(sig2, m)
    if(length(sig2) != m || !is.numeric(sig2)|| any(sig2 <=0))
	stop("'sig2' = sigma^2	must be > 0 with same length as 'mu'")
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
    structure(name = name, class = "norMix",
	      .Data = cbind(mu = mu, sig2 = sig2, w = w))
}

is.norMix <- function(obj)
{
  ## Purpose: is 'obj' a "norMix", i.e.	 Normal Mixture object ?
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:38
  inherits(obj, "norMix") &&
  (!is.null(w <- obj[,"w"])) &&
  is.numeric(w) && all(w >= 0) && abs(sum(w)-1) < 1000*.Machine$double.eps
}

m.norMix <- function(obj) nrow(obj) ##	Number of components of	 normal mixture

mean.norMix <- function(x, ...)
{
  ## Purpose: Return "true mean", i.e., the expectation of  a normal mixture.
  if(!is.norMix(x)) stop("'x' must be a 'Normal Mixture' object!")
  x[,"w"] %*% x[,"mu"]
}

var.norMix <- function(x, ...)
{
  ## Purpose: 'true' Variance, i.e. E[(X- E[X])^2]  for X ~ normal mixture.
  if(!is.norMix(x)) stop("'x' must be a 'Normal Mixture' object!")
  w <- x[,"w"]
  mj <- x[,"mu"]
  mu <- w %*% mj
  w %*% (x[,"sig2"] + (mj - mu)^2)
}

print.norMix <- function(x, ...)
{
    ## Purpose: print method for  "norMix" objects (normal mixtures)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 97, 10:02
    ox <- x
    has.nam <- !is.null(nam <- attr(x,"name"))
    cat("'Normal Mixture' object",
	if(has.nam) paste("\t ``", nam, "''", sep=''), "\n")
    if(has.nam) attr(x, "name") <- NULL
    cl <- class(x);  cl <- cl[ cl != "norMix"] #- the remaining classes
    class(x) <- if(length(cl) > 0) cl ## else NULL
    NextMethod("print", ...)
    invisible(ox)
}

sort.norMix <- function(x, decreasing = FALSE, ...) {
    ## sort according to 'mu'
    x[] <- x[sort.list(x[,"mu"], decreasing = decreasing, ...) , ]
    x
}

dnorMix <- function(x, obj, log = FALSE)
{
  ## Purpose: density evaluation for "norMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Arguments: x: numeric; obj: Normal Mixture object
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj)) {
      ## Old version had (obj, x, ..):
      if(is.norMix(x)) { ## swap the first two arguments
          tmp <- x ; x <- obj; obj <- tmp
          warning("Deprecated use of dnorMix(obj, x, ..);
  Either use dnorMixL(), or the new argument order (x, obj, ...) and
  note that dnorMix() returns a numeric vector (not a list).")
      }
      else stop("'obj' must be a 'Normal Mixture' object!")
  }
  w <- obj[,"w"]; mu <- obj[,"mu"]; sd <- sqrt(obj[,"sig2"])
  m <- length(w) #-- number of components
  if(m == 1)
      return(dnorm(x, mean = mu[1], sd = sd[1], log = log))
  ## else
  y <- numeric(length(x))
  for(j in 1:m)
    y <- y + w[j] * dnorm(x, mean = mu[j], sd = sd[j])
  if(log) log(y) else y
}

dnorMixL <- function(obj, x = NULL, log = FALSE, xlim = NULL, n = 511)
{
  ## Purpose: density evaluation for "norMix" objects (normal mixtures)
  ## -------------------------------------------------------------------------
  ## Arguments: obj: Normal Mixture object;  x:
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj))
    stop("'obj' must be a 'Normal Mixture' object!")
  if(is.null(x)) {
    if(is.null(xlim)) ##-- construct "reasonable" abscissa values
      xlim <- mean.norMix(obj) + c(-3,3)*sqrt(var.norMix(obj))
    x <- seq(xlim[1], xlim[2], length = n)
  }
  list(x = x, y = dnorMix(x, obj, log=log))
}


rnorMix <- function(n, obj)
{
    ## Purpose: Generate random numbers according to "norMix" object `obj'
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Jun 2002, 16:03
    mu <- obj[,"mu"]
    sd <- sqrt(obj[,"sig2"])
    if(n == 1) {
	j <- sample(length(mu), size = 1, prob = obj[,"w"])
	rnorm(1, mean = mu[j], sd = sd[j])
    } else {
	nj <- as.vector(rmultinom(n=1, size = n, prob = obj[,"w"]))
	sample(unlist(lapply(seq(along=nj), function(j)
			     rnorm(nj[j], mean = mu[j], sd = sd[j]))))
    }
}

## From: Erik Jørgensen <Erik.Jorgensen@agrsci.dk>
## Date: Thu, 13 Nov 2003 02:06:27 +0100
##
## ....... Please, feel free to use them.
##
## Erik Jørgensen
## Danish Institute of Agricultural Sciences

pnorMix <- function(q, obj, lower.tail = TRUE, log.p = FALSE)
{
    if (!is.norMix(obj)) {
        ## Old version had (obj, q):
        if(is.norMix(q)) { ## swap the first two arguments
            tmp <- q ; q <- obj; obj <- tmp
            warning("Deprecated use of pnorMix(obj, q, ..); NEW argument order is (q, obj, ...)")
        }
        else stop("'obj' must be a 'Normal Mixture' object!")
    }
    sd <- sqrt(obj[,"sig2"])
    ## if log.p just log(.) at the end [to be more accurate, need much more..]
    cc <- if(log.p) function(m) log(c(m)) else c
    ## q can be a vector: -> outer
    cc(pnorm(sweep(outer(q, obj[,"mu"], "-"), 2, sd, "/"),
             lower.tail= lower.tail) %*% obj[, "w"])
}

dpnorMix <- function(x, obj, lower.tail = TRUE)
{
    ## Purpose: compute dnorMix() and pnorMix() simultaneously
    ## ----------------------------------------------------------------------
    ## Arguments: x: numeric; obj: 'norMix' object
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 3 Jan 2008
    stopifnot (is.norMix(obj))
    w <- obj[,"w"]; mu <- obj[,"mu"]; sd <- sqrt(obj[,"sig2"])
    ## This looks smarter, but really is slower :
    ##     z <- sweep(outer(x, obj[,"mu"], "-"), 2, sqrt(obj[,"sig2"]), "/")
    ##     list(d = c(dnorm(z) %*% w),
    ##          p = c(pnorm(z, lower.tail= lower.tail) %*% w))
    m <- length(w) #-- number of components
    d <- p <- numeric(length(x))
    for(j in 1:m) {
        d <- d + w[j] * dnorm(x, mean= mu[j], sd= sd[j])
        p <- p + w[j] * pnorm(x, mean= mu[j], sd= sd[j], lower.tail=lower.tail)
    }
    list(d = d, p = p)
}

qnorMix <-
    function(p, obj, lower.tail = TRUE, log.p = FALSE,
	     tol = .Machine$double.eps^0.25, maxiter = 1000, traceRootsearch = 0,
	     method = c("interpQspline", "interpspline", "eachRoot", "root2"),
	     l.interp=20)
    ## NOTE: keep defaults consistent with 'uniroot':
{
  if (!is.norMix(obj)) {
    ## Old version had (obj, p):
    if(is.norMix(p)) { ## swap the first two arguments
      tmp <- p ; p <- obj; obj <- tmp
      warning("Deprecated use of qnorMix(obj, p, ..); NEW argument order is (p, obj, ...)")
    }
    else stop("'obj' must be a 'Normal Mixture' object!")
  }
  mu <- obj[, "mu"]
  sd <- sqrt(obj[, "sig2"])
  m <- m.norMix(obj)
  if(m == 1) # one component
      return(qnorm(p, mu, sd, lower.tail=lower.tail, log.p=log.p))

  ## else

  S <- if(lower.tail) 1 else -1
  ## vectorize in 'p'
  r <- p
  ## Solve the left/right extremes p \in {0 , 1}
  left <- if(log.p) rep(FALSE, length(p)) else p <= 0
  right <- p >= if(log.p) 0 else 1
  r[left] <- -Inf*S ; r[right] <- Inf*S
  imid <- which(mid <- !left & !right) # 0 < p < 1
  if(length(imid)) {

      f.make <- function(p.i) {
	  if(traceRootsearch >= 3)
	      function(l) {
		  p <- pnorMix(l, obj,
			       lower.tail=lower.tail, log.p=log.p)
		  cat(sprintf("p(%-19.16g) = %-19.16g\n", l, p))
	      }
	  else
	      function(l) pnorMix(l, obj, lower.tail=lower.tail,
				  log.p=log.p) - p.i
      }
      outRange <- function(p.i)
	  range(qnorm(p.i, mu, sd, lower.tail=lower.tail, log.p=log.p))

      safeUroot <- function (f, interval,
			     lower = min(interval), upper = max(interval),
			     tol=tol, maxiter=maxiter, ...)
      {
	  if(traceRootsearch >= 2)
	      cat(sprintf("search in [%g,%g]\n", lower, upper))

	  ## make sure we have S*f(lower) < 0 and S*f(upper) > 0:
	  delta.r <- 0.01*max(1e-7, abs(lower))
	  while(S*(f.lo <- f(lower)) > 0) {
	      lower <- lower - delta.r
	      if(traceRootsearch)
		  cat(sprintf(" .. modified lower: %g\n",lower))
	      delta.r <- 2 * delta.r
	  }
	  delta.r <- 0.01*max(1e-7, abs(upper))
	  while(S*(f.up <- f(upper)) < 0) {
	      upper <- upper + delta.r
	      if(traceRootsearch)
		  cat(sprintf(" .. modified upper: %g\n",upper))
	      delta.r <- 2 * delta.r
	  }

	  ## Here, we require R >= 2.6.0 with the improved uniroot():
	  uniroot(f, lower=lower, upper=upper,
		  f.lower = f.lo, f.upper = f.up,
		  tol=tol, maxiter=maxiter, ...)
      }

      ## sort p[] increasingly for easier root finding start:
      p <- sort(p[mid], index.return = TRUE)
      ip <- imid[p$ix]
      pp <- p$x

      hasDup <- any(iDup <- duplicated(pp))
      if(hasDup) {
	  isUniq <- !iDup
	  ## want *strictly* increasing sometimes; save CPU anyway
	  pp <- pp[isUniq]
	  i.pp <- cumsum(isUniq) ## pp[i.pp]  |-->  original pp[]
      }

      np <- length(pp)
      rr <- pp #- rr will contain = q..mix(pp, *)
      method <- if(np <= 2) "eachRoot" else match.arg(method)

      if(method == "eachRoot") { ## root finding from left to right ...
	  for(i in seq(along=pp)) {
	      ff <- f.make(pp[i])
	      rq <- outRange(pp[i])
	      ## since pp[] is increasing, we can start from last 'root':
	      if(i > 1 && rq[1] < root)
		  rq[1] <- root
	      root <- safeUroot(ff, interval = rq, tol=tol, maxiter=maxiter)$root
	      rr[i] <- root
	  }
      }
      else { ## => np > 2
	rr[1] <- safeUroot(f.make(pp[1]), interval = outRange(pp[1]),
                           tol=tol, maxiter=maxiter)$root
        rr[np] <- safeUroot(f.make(pp[np]), interval = outRange(pp[np]),
                            tol=tol, maxiter=maxiter)$root
        ni <- length(iDone <- as.integer(c(1,np)))
        if(any(method == c("interpQspline", "interpspline"))) {
            ## reverse interpolate, using relatively fast pnorMix()!

            pp. <- pp[-iDone]
            ## those mu's that are inside our range:
            rXtr <- rr[if(lower.tail) c(1L,np) else c(np,1L)]
            mu. <- unique(sort(mu[rXtr[1] < mu & mu < rXtr[2]]))
            k <- length(qs <- c(rXtr[1], mu., rXtr[2]))
            qs. <- qs[-k]
            dq <- qs[-1] - qs.      # == delta(qs)
            ## l.interp values between each mu
            stopifnot(l.interp >= 1)
            qi <- c(t(dq %*% t(seq_len(l.interp)/l.interp) + qs.))
            stopifnot(!is.unsorted(qi)) ## << FIXME remove if never triggering
            ppi <- pnorMix(qi, obj, lower.tail=lower.tail, log.p=log.p)

            ## in an extreme case, pnorMix() is horizontal; hence
            ## qnorMix() has practically a discontinuity there.
            ## In that case, splinefun() completely "fails";
            ## we do need  *monotone* (spline) interpolation:
            mySfun <- {
                if(getRversion() < "2.8.0") monoHsplinefun
                else function(x, y, ...) splinefun(x,y, ..., method="monoH.FC")
            }
            if(method == "interpspline")
                qpp <- mySfun(ppi, qi)(pp.) ## is very fast
            else {                          ## "interpQspline
                ## logit() transform the P's --> interpolation is more linear
                muT <- c(obj[, "w"] %*% mu)
                qpp <- mySfun(qlogis(ppi, muT), qi)(qlogis(pp., muT))
            }

            if(log.p)
                warning("Newton steps for 'log.p = TRUE' not yet implemented")
            else {
                ## now end with a few Newton steps
                for(k in 1:maxiter) {
                    dp <- dpnorMix(qpp, obj, lower.tail=lower.tail)
                    del.p <- dp$p - pp.
                    ## FIXME?: del.p  may suffer from considerable cancellation
                    relE.f <- abs(del.p)
                    n0 <- relE.f > 0 & pp. > 0
                    relE.f[n0] <- (relE.f/pp.)[n0]
                    ii. <- dp$d > 0 ## & relE.f > tol
                    if(!any(ii.)) {
                        relErr <- mean(relE.f)
                        break           # not converged though
                    }
                    ## del.q := Delta(q) =  F(q) / f(q)   or 0 if f(q)=0
                    del.q <- numeric(length(pp.))
                    del.q[ii.] <- S*(del.p/dp$d)[ii.]
                    ## only modify qpp[] where Newton step is ok:
                    ## e.g. resulting qpp must remain increasing
                    while(length(iF <- which(S*diff(qNew <- qpp - del.q) <= 0))) {
                        iF <- c(iF,iF+1L)
                        del.q[iF] <- del.q[iF] / 2
                        if(traceRootsearch) cat(",")
                    }
                    qpp[ii.] <- qNew[ii.]
                    relErr <- sum(abs(del.q[ii.])) / sum(abs(qpp[ii.]))
                    if(traceRootsearch) {
                        cat(k,": relE =", formatC(relErr), sep='')
                        if(traceRootsearch >= 2) {
                            cat(" |\n")
                            if(traceRootsearch == 2 || length(qpp) <= 10)
                                print(summary(del.q / abs(qpp)))
                            else print(del.q / abs(qpp))
                        }
                        else cat("\n")
                    }
                    if(relErr < tol) break
                }
                if(relErr >= tol)
                    warning("Newton iterations have not converged")
            }
            rr[-iDone] <- qpp
        }
	  ## method == "root2"
	  else while(ni < np) {
	      ## not "done";  ni == length(iDone)
	      oi <- iDone
	      i.1 <- oi[-ni]
	      i.2 <- oi[-1]
	      l.new <- i.2 > i.1 + 1L	# those "need new"
	      ii <- which(l.new)
	      iN <- (i.1 + i.2 + 1L) %/% 2
	      stopifnot( (i.1 < iN)[ii], (iN < i.2)[ii])
	      if(traceRootsearch) cat("ni new intervals, ni=",ni,"\n")
	      for(j in ii) {
		  ## look in between i.1[j] .. i.2[j]
		  ## NB: we can prove that  i.1[j] < iN[j] < i.2[j]
		  rr[iN[j]] <- safeUroot(f.make(pp[iN [j]]),
					 lower= rr[i.1[j]],
					 upper= rr[i.2[j]],
					 tol=tol, maxiter=maxiter)$root
	      }
	      ## update iDone[]:
	      seq_old <- seq_len(ni)
	      ni <- ni + length(ii)
	      iDone <- integer(ni)
	      iDone[seq_old	  + c(0L, cumsum(l.new))] <- oi
	      iDone[seq_along(ii) + ii			] <- iN[ii]
	  }
      } # else { method ..}

      r[ip] <- if(hasDup) rr[i.pp] else rr

  } # end if(np)
  r
}## end{qnorMix}

plot.norMix <-
    function(x, type = "l", n = 511, xout = NULL, xlim = NULL,
	     xlab = "x", ylab = "f(x)", main = attr(x,"name"), lwd = 1.4,
	     p.norm = TRUE, p.h0 = TRUE, p.comp = FALSE,
	     parNorm = list(col= 2, lty = 2, lwd = 0.4),
	     parH0   = list(col= 3, lty = 3, lwd = 0.4),
	     parComp = list(col= "blue3", lty = 3, lwd = 0.4),
	     ...)
{
    ## Purpose: plot method for	 "norMix" objects (normal mixtures)
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 20 Mar 1997
    if(!is.null(xlim) && is.null(xout)) ## determine xout
	xout <- seq(xlim[1], xlim[2], length = n)
    d.o <- dnorMixL(x, x = xout, n = n)
    if(p.norm)
	dn <- dnorm(d.o$x, mean = mean.norMix(x), sd = sqrt(var.norMix(x)))
    if(!is.null(ll <- list(...)[["log"]]) && "y" %in% strsplit(ll,"")[[1]])
	y0 <- max(1e-50, min(d.o$y, if(p.norm) dn))
    else y0 <- 0
    plot(d.o, type = type, xlim = xlim, ylim = c(y0, max(d.o$y, if(p.norm) dn)),
	 main = main, xlab = xlab, ylab = ylab, lwd = lwd, ...)
    if(p.norm)	do.call("lines",  c(list(x = d.o$x, y = dn), parNorm))
    if(p.h0)	do.call("abline", c(list(h = 0), parH0))
    if(p.comp) {
        m <- m.norMix(x) #-- number of components
        w <- x[,"w"]; mu <- x[,"mu"]; sd <- sqrt(x[,"sig2"])
        for(j in 1:m)
            do.call("lines",
                    c(list(x = d.o$x,
                           y = w[j] * dnorm(d.o$x, mean = mu[j], sd = sd[j])),
                      parComp))
    }
    invisible(x)
}

lines.norMix <-
    function(x, type = "l", n = 511, xout = NULL, lwd = 1.4,
	     p.norm = FALSE, parNorm = list(col = 2, lty = 2, lwd = 0.4), ...)
{
    ## Purpose: lines method for "norMix" objects (normal mixtures)
    ## -------------------------------------------------------------
    ## Author: Martin Maechler, Date: 27 Jun 2002, 16:10
    xlim <- if(is.null(xout)) par("usr")[1:2] # else NULL
    d.o <- dnorMixL(x, x = xout, n = n, xlim = xlim)
    lines(d.o, type = type, lwd = lwd, ...)
    if(p.norm) {
	dn <- dnorm(d.o$x, mean = mean.norMix(x), sd = sqrt(var.norMix(x)))
	do.call("lines", c(list(x = d.o$x, y = dn), parNorm))
    }
    invisible()
}


r.norMix <- function(obj, x = NULL, xlim = NULL, n = 511, xy.return = TRUE)
{
  ## Purpose: Compute r := f / f0; f = normal mixture; f0 = "best" normal approx
  ## Author : Martin Maechler, Date: 20 Mar 97, 10:14
  if(!is.norMix(obj)) stop("'obj' must be a 'Normal Mixture' object!")
  d.o <- dnorMixL(obj, x, xlim = xlim, n = n)
  dn  <- dnorm(d.o$x, mean = mean.norMix(obj), sd = sqrt(var.norMix(obj)))
  if(xy.return) list(x = d.o$x, y = d.o$y / dn, f0 = dn) else d.o$y / dn
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


nM2par <- function(obj)
{
    ## Purpose: translate norMix object into our parametrization par.vector
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 17 Dec 2007
    stopifnot(is.norMix(obj))
    ## logit() == qlogis(); log(sqrt(.)) = log(.)/2
    c(qlogis(obj[-1,"w"]), obj[,"mu"], log(obj[,"sig2"])/2)
}
## FIXME? This is nowhere used
.nM2par <- function(lst)
{
    ## Purpose: Fast version of nM2par()
    ## -------------------------------------------------
    ## Author: Martin Maechler, Date: 18 Dec 2007
    c(qlogis(lst$w[-1]), lst$mu, log(lst$sig2)/2)
}


par2norMix <- function(p, name = sprintf("{from %s}",
			  deparse(substitute(p))[1]))
{
    ## Purpose: build norMix object from our parametrization par.vector
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 17 Dec 2007

    force(name) # substitute(..)
    lp <- length(p)
    stopifnot(is.numeric(p), lp %% 3 == 2)
    m <- (lp + 1L) %/% 3
    m1 <- m - 1L
    names(p) <- NULL # so they are not transferred to mu,...
    mu  <- p[m:(m+m1)]
    sig2 <- exp(2 * p[(m+m):(m+m+m1)]) ## sigma^2 = exp(tau)^2 = exp(2*tau)
    if(m == 1)
        norMix(mu, sig2, name= name)
    else { ## -- m >= 2
        pi. <- plogis(p[1:m1]) ## \pi_j = inv_logit(\lambda_j)
        if((sp <- sum(pi.)) > 1)
            stop(sprintf("weights sum up to %.3g > 1 !", sp))

        norMix(mu, sig2, w = c(1 - sp, pi.), name = name)
    }
}


if(FALSE) ## this is not needed -- but mention it on ?llnorMix
logLiknorMix <- function(obj, x) {
    ## Purpose: log-likelihood for 'norMix'
    sum(dnorMix(x, obj, log=TRUE))
}

llnorMix <- function(p, x, m = (length(p)+1)/3)
{
    ## Purpose: log-likelihood
    ## ----------------------------------------------------------------------
    ## Arguments: p : parameter vector, see below
    ##            x : data vector
    ##            m : number of mixture components
    ##
    ##  'p' is particularly parametrized:
    ##		  p = c( lambda_j, mu_j, tau_j)	 where
    ##		\lambda_j = logit(\pi_j), j=2,..,m; and \pi_1 := 1- sum_j\pi_j
    ##	    and \tau_j = log(\sigma_j)	such that parameters are unconstrained
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 17 Dec 2007, 17:41
    stopifnot(is.numeric(x), is.numeric(p), !is.matrix(p), 3*m == length(p)+1,
              m == (m. <- as.integer(m)), (m <- m.) >= 1)
    m1 <- m - 1L
    mu	<- p[m:(m+m1)]
    sigma <- exp(p[(m+m):(m+m+m1)]) ## sigma = exp(tau)
    if(m == 1)
	return( sum(dnorm(x, mean = mu[1], sd = sigma[1], log = TRUE)) )

    ## else -- m >= 2
    pi. <- plogis(p[1:m1]) ## \pi_j = inv_logit(\lambda_j)
    if((sp <- sum(pi.)) > 1) return(-Inf) # worst possible value
    ## pi. <- c(pi., 1 - sp) # \pi_1 := 1 - sum_{j=2}^{m} \pi_j
    y <- (1 - sp) * dnorm(x, mean = mu[1], sd = sigma[1])
    for(j in 2:m)
	y <- y + pi.[j- 1L] * dnorm(x, mean = mu[j], sd = sigma[j])
    ## return
    sum(log(y))
}


clus2norMix <- function(gr, x, name = deparse(sys.call()))
{
    ## Purpose: "Clustering to normal Mixture"
    ## ----------------------------------------------------------------------
    ## Arguments: gr: grouping/clustering vector (in {1,..,K}); possibly factor
    ##		  x : (original) data vector
    ##		name : name for norMix() object; by default constructed
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 31 Dec 2007, 10:24

    if(length(gr) != (n <- length(x)))
	stop("'gr' and 'x' are not of the same length")
    r <- lapply(unname(split(x,gr)), ## << simple version of tapply(x, gr, *)
		function(u){ n <- length(u); m <- mean(u)
			     list(m, sum((u - m)^2)/(n-1), n) })
    norMix(mu	= unlist(lapply(r,`[[`, 1L)),
	   sig2 = unlist(lapply(r,`[[`, 2L)),
	   w	= unlist(lapply(r,`[[`, 3L))/n,
	   name = name)
}
