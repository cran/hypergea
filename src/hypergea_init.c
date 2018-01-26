#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void hypergeom_2x2x2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hypergeom_IxJ(void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP lfactorial(SEXP);

static const R_CMethodDef CEntries[] = {
    {"hypergeom_2x2x2", (DL_FUNC) &hypergeom_2x2x2, 9},
    {"hypergeom_IxJ",   (DL_FUNC) &hypergeom_IxJ,   9},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
 //   {"hypergeom_IxJ_a", (DL_FUNC) &hypergeom_IxJ_a, 6},
    {"lfactorial",      (DL_FUNC) &lfactorial,      1},
    {NULL, NULL, 0}
};

void R_init_hypergea(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
