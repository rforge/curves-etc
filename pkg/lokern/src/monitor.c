// Flexible printing of informative messages from Fortran
// -----------------                         ------------
//-- original: ~/R/Pkgs/cobs99/src/monitor.c

#include <R.h>

/* called for  trace >= 1 : ----------------------------------------------- */

void F77_SUB(monit0)(int *kind, int *n, int *m, int *nue, int *kord,
		     int *inputb, int *isrand, double *b, int *trace)
{
    if(*inputb)
	Rprintf("%s (n=%d,m=%d; nue=%d, isrand=%d, inputp=%d, b=%g) -> kord=%d\n",
		(*kind == 0)? "glkerns" : "lokerns",
		*n, *m, *nue, *isrand, *inputb, *b, *kord);
    else
	Rprintf("%s (n=%d,m=%d; nue=%d, isrand=%d, inputp=%d) -> kord=%d\n",
		(*kind == 0)? "glkerns" : "lokerns",
		*n, *m, *nue, *isrand, *inputb,     *kord);
}

void F77_SUB(monit1)(int *phase, int *trace)
{
    Rprintf(" phase %2d\n", *phase);
    // TODO print more when (*trace >= 2) ...
}

// from the end of "loop 100" of glkerns():
void F77_SUB(monit_s)(double *r2, double *osig, double *q, double *sig) {
    Rprintf("  r2 = %g, old sig = %g, q = sig/rvar = %g, new sig = %g --> goto 100\n",
	    *r2, *osig, *q, *sig);
}



// Called only from kernel() and kernp()  in auxkerns.f
void F77_SUB(monitk0)(int *caller, int *n, int *m, double *b, double *chan,
		      double *chR, int *do_classic)
{
    Rprintf(" %s(n=%3d,m=%3d; b=%9.7g) -> (chg.pt,cut_b)=(%4.1f,%5.2f) => '%s'\n",
	    // caller = 0 | 1  <==> from (kernel | kernp )
	    (*caller) ? "kernp" : "kernel",
	    *n, *m, *b, *chan, *chR, (*do_classic) ? "classic" : "fast O(n)");
}

// Called from kernfa() and kernfp() in auxkerns.f
void F77_SUB(monitfp)(int *caller, int *n, double *b, int *nue,
		      int *kord, int *ny, int *m, int *trace)
{
    Rprintf("  %s(n=%3d, b=%9.7g, nue=%d, kord=%d, ny=%d, m=%d):\n",
	    // caller = 0 | 1  <==> from ( kernfa | kernfp )
	    (*caller) ? "kernfp" : "kernfa",
	    *n, *b, *nue, *kord, *ny, *m);
}


/* called for  trace >= 2 : ----------------------------------------------- */

// inside smo() :
void F77_SUB(monits)(double *tau, int *ist, int *n, int *iboun)
{
    Rprintf("   smo(t=%12.6g, j = ist:n = %d:%d); iboun=%2d\n",
	    *tau, *ist, *n, *iboun);
}
