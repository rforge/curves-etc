#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
// glkerns.f :
void F77_NAME(glkern_s)(double *t, double *x, double *tt, double *y, int *n, int *m, int *nue, int *kord, int *hetero, int *isrand, int *inputb, int *m1, double *tl, double *tu, double *s, double *sig, double *wn, double *w1, double *b, int *trace);
// lokerns.f :
void F77_NAME(lokern_s)(double *t, double *x, double *tt, double *y, int *n, int *m, int *nue, int *kord, int *hetero, int *isrand, int *inputb, int *m1, double *tl, double *tu, double *s, double *sig, double *wn, double *w1, double *wm, double *ban, int *trace);
// auxkerns.f :
void F77_NAME(resest)(double *t, double *x, int *n, double *res, double *snr, double *sigma2);

static const R_FortranMethodDef FortranEntries[] = {
    {"glkern_s", (DL_FUNC) &F77_SUB(glkern_s), 20},
    {"lokern_s", (DL_FUNC) &F77_SUB(lokern_s), 21},
    {"resest",   (DL_FUNC) &F77_SUB(resest),    6},
    {NULL, NULL, 0}
};

void R_init_lokern(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
