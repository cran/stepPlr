#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(solveplr)(double *x, int *n, double *z, int *lenz, int *nobs, double *lam, int *status);

static const R_FortranMethodDef FortranEntries[] = {
    {"solveplr", (DL_FUNC) &F77_NAME(solveplr), 7},
    {NULL, NULL, 0}
};

void R_init_stepPlr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
