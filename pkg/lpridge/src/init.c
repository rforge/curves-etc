#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "lpridge.h"
// -> lpepa_s , lpridge_s

static const R_FortranMethodDef FortranEntries[] = {
    {"lpepa_s",   (DL_FUNC) &F77_NAME(lpepa_s),   16},
    {"lpridge_s", (DL_FUNC) &F77_NAME(lpridge_s), 20},
    {NULL, NULL, 0}
};

void R_init_lpridge(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
