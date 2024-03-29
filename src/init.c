#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _bcROCsurface_asy_var_vus_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _bcROCsurface_vus_c(SEXP, SEXP);
extern SEXP _bcROCsurface_vus_c_full(SEXP, SEXP, SEXP);
extern SEXP _bcROCsurface_vus_c_full_core(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_bcROCsurface_asy_var_vus_c",  (DL_FUNC) &_bcROCsurface_asy_var_vus_c,   8},
  {"_bcROCsurface_vus_c",          (DL_FUNC) &_bcROCsurface_vus_c,           2},
  {"_bcROCsurface_vus_c_full",     (DL_FUNC) &_bcROCsurface_vus_c_full,      3},
  {"_bcROCsurface_vus_c_full_core",(DL_FUNC) &_bcROCsurface_vus_c_full_core, 3},
  {NULL, NULL, 0}
};

void R_init_bcROCsurface(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
