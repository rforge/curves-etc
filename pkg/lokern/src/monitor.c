// Flexible printing of informative messages from Fortran
// -----------------                         ------------
//-- original: ~/R/Pkgs/cobs99/src/monitor.c

#include <R.h>

/* called for  trace >= 1 : ----------------------------------------------- */

void F77_SUB(monit0)(int *kind, int *n, int *m, int *nue, int *kord, int *trace)
{
    Rprintf("%s (n=%d,m=%d; nue=%d) -> kord=%d\n",
	    (*kind == 0)? "glkerns" : "lokerns",
	    *n,*m, *nue, *kord);
    // TODO print more when (*trace >= 2) ...
}

void F77_SUB(monit1)(int *phase, int *trace)
{
    Rprintf(" phase %2d\n", *phase);
    // TODO print more when (*trace >= 2) ...
}

void F77_SUB(monitk0)(int *n, int *m, double *b, double *chan,
		      double *chR, int *do_classic)
{
    Rprintf(" kernel(n=%3d,m=%3d; b=%12.7g) -> (chg.pt,cut_b)=(%4.1f,%5.2f) => '%s'\n",
	    *n, *m, *b, *chan, *chR, (*do_classic) ? "classic" : "fast O(n)");
}


/* called for  trace >= 2 : ----------------------------------------------- */

// inside smo() :
void F77_SUB(monits)(double *tau, int *ist, int *n, int *iboun)
{
    Rprintf("   smo(t=%12.6g, j = ist:n = %d:%d); iboun=%2d\n",
	    *tau, *ist, *n, *iboun);
}
