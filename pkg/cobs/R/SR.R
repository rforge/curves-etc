### --> ../convexreg-PietGr/  for the original code, and more

## FIXME: copied the arg list (and start) from conreg().
## -----  BUT  SR() does not (yet) work with weights
SR <- function(x, y = NULL, w = NULL, convex = FALSE,
	       tol = c(1e-10, 1e-10), maxit = c(1000, 20),
               ## conreg has (200, 20) -- to check/test
	       adjTol = TRUE, verbose = FALSE)
{
  ## Author: Martin Maechler, Date: 27 Jul 2012, 21:35

  ## --- use the same code as in smooth.spline() {!}
  xy <- xy.coords(x, y)
  y <- xy$y
  x <- xy$x
  if(!is.null(w))
      stop("weights are not yet supported.  Use conreg() for that.")
  ## n <- length(x)

  ## cobs() has all the weights; handling of ties, etc etc
  ## Here, we only sort (if needed)
  if(doSort <- is.unsorted(x)) {
      i. <- sort.list(x, method="quick")
      x. <- x[i.]
      y. <- y[i.]
  } else {
      x. <- x
      y. <- y
  }

  n <- length(x.)
  x.rng <- x.[n] - x.[1L]

  if(!convex) y <- - y # and will revert at the end

  stopifnot(length(maxit) >= 1, maxit == round(maxit))
  if(length(maxit) == 1) maxit <- rep.int(maxit, 2)
  stopifnot(length(tol) >= 1, tol >= 0)
  if(length(tol) == 1) tol <- rep.int(tol, 2)

  ## rescale 'tol' to be x-scale equivariant:
  if(n > 1) tol <- tol / sd(x.)

  r <- .C("SR_R",
	  n=  as.integer(n),
	  cc= as.double(x.rng),
	  m1= integer(1),
	  ind= integer(n+1),
	  x = as.double(x),
	  y = as.double(y),
	  r = double(n),
	  ## MM: what are these? -- surely some are interesting
	  ##  (if not, create and free in C !!)
	  R = double(n),
	  H = double(n),
	  S = double(n),
	  Y = double(n),
	  D = double(n),
	  tol = as.double(tol),
	  maxit = as.integer(maxit),
	  phiBL = double(1),
	  numIt = integer(1),
	  PACKAGE="cobs")
	## the vector r contains the convex LS estimates f[i] at the points x[i]
	## the data are in the vector y and phi is the criterion function
	## the vector of indices ind contains the indices of the points where the
	## solution has a kink (only for method 1) and m is the number of kinks

  if(convex) {
      M <- r$ r
      H <- r$ H
  } else { ##  (!convex)	switch the signs back
      y. <- -y.
      M <- -r$ r
      H <- -r$ H
      ## Conv <- -Conv
  }

    if(doSort) {
        ## FIXME -- revert or ?????
    }

  ## return
  structure(list(x = x., y = y., ##  w = w.,
		 yf = M,
		 convex = convex, call = match.call(),
		 iKnots = r$ind[1L + seq_len(r$m1)],
		 deriv.loc = H,
		 ##FIXME conv.loc = Conv,
		 iter = r$numIt, # c(iter, innerIt)
		 phi = r$phiBL
		 ),
	    class = "conreg") ## <<- find a better name
}
