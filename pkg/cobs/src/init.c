#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "splines.h"
/* other .C calls */
/* extern void SR_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *); */

static const R_CMethodDef CEntries[] = {
    {"spline_basis", (DL_FUNC) &spline_basis, 8},
    {"spline_value", (DL_FUNC) &spline_value, 8},
    /* {"SR_R", (DL_FUNC) &SR_R, 16}, */
    {NULL, NULL, 0}
};

void R_init_cobs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
