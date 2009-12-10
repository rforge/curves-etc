###  safeUroot()  used to be  *internal*  to qnorMix()
###  -----------  but as it's generally useful,
##   {pay a (very?) small penalty and make it public __ TODO: currently NAMESPACE-hidden!
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
