// Flexible printing of informative messages from Fortran
// -----------------                         ------------
//-- original: ~/R/Pkgs/cobs99/src/monitor.c

#include <R.h>

/* called for  trace >= 1 : ----------------------------------------------- */

void F77_NAME(monit0)(int *kind, int *n, int *m, int *nue, int *kord,
		     int *inputb, int *isrand, double *b, int *trace)
{
    if(*inputb)
	Rprintf("%s (n=%d,m=%d; nue=%d, isrand=%d, inputp=%d, b=%10.7g) -> kord=%d\n",
		(*kind == 0)? "glkerns" : "lokerns",
		*n, *m, *nue, *isrand, *inputb, *b, *kord);
    else
	Rprintf("%s (n=%d,m=%d; nue=%d, isrand=%d, inputp=%d) -> kord=%d\n",
		(*kind == 0)? "glkerns" : "lokerns",
		*n, *m, *nue, *isrand, *inputb,     *kord);
}

void F77_NAME(monit1)(int *phase, int *trace)
{
    Rprintf(" phase %2d\n", *phase);
    // TODO print more when (*trace >= 2) ...
}

// after label '40' in phase 4 of glkerns() / lokerns()
void F77_NAME(monit40)(int* itt, int*  il, int* iu, double *tl, double *tu, int *trace)
{
    Rprintf(" 40: itt=%d, (il=%d, iu=%d) -  (tl=%g, tu=%g)\n",
	    *itt, *il,*iu, *tl,*tu);
    // TODO print more when (*trace >= 2) ...
}

// after label '100' in phase 10 of {gl,lo}kerns() : monit1c(il,iu, vi, trace)
void F77_NAME(monit1c)(int* il, int* iu, double *vi, int *trace)
{
    Rprintf("  100: (il=%d, iu=%d) -> vi=%g\n", *il,*iu, *vi);
    // if(*trace >= 2) <print more>
}

// after label '111' in phase 11 of {gl,lo}kerns():   monit111(i, s(i-1),ssi)
void F77_NAME(monit111)(int* i, double *s1, double *s2)
{
    Rprintf("  111: i=%d: swapped s[i-1] w/ s[i]: (%g, %g)\n", *i, *s1, *s2);
}


// from the end of "loop 100" of glkerns():
void F77_NAME(monit_s)(double *r2, double *osig, double *q, double *sig) {
    Rprintf("  r2 = %g, old sig = %g, q = sig/rvar = %g, new sig = %g --> goto 100\n",
	    *r2, *osig, *q, *sig);
}

void F77_NAME(monit2ib)(int *trace) {
    Rprintf(" (inputb && !isrand) --> goto 123\n");
}


// Called only from kernel() and kernp()  in auxkerns.f
void F77_NAME(monitk0)(int *caller, int *n, int *m, double *b, int *ny,
		       double *chan, double *chR, int *do_classic)
{
    Rprintf(" %s(n=%3d, m=%3d, ny=%d => %s; b=%10.7g) -> (chg.pt,cut_b)=(%4.1f,%5.2f) => '%s'\n",
	    // caller = 0 | 1  <==> from (kernel | kernp )
	    (*caller) ? "kernp" : "kernel",
	    *n, *m, *ny, (*ny) ? "local bandw." : "global bandw.",
	    *b, *chan, *chR, (*do_classic) ? "classic" : "fast O(n)");
}

// Called from kernfa() and kernfp() in auxkerns.f
void F77_NAME(monitfp)(int *caller, int *n, double *b, int *nue,
		       int *kord, int *ny, int *m, int *trace)
{
    Rprintf("  %s(n=%3d, b=%10.7g, nue=%d, kord=%d, ny=%d, m=%d):\n",
	    // caller = 0 | 1  <==> from ( kernfa | kernfp )
	    (*caller) ? "kernfp" : "kernfa",
	    *n, *b, *nue, *kord, *ny, *m);
}


/* called for  trace >= 2 : ----------------------------------------------- */

// inside smo() :
void F77_NAME(monits)(double *tau, int *ist, int *n, int *iboun)
{
    Rprintf("   smo(t=%12.6g, j = ist:n = %d:%d); iboun=%2d\n",
	    *tau, *ist, *n, *iboun);
}
