#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "mgcv.h"

extern void gsl_bspline(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gsl_bspline_deriv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RuniqueCombs(void *, void *, void *, void *);

R_CallMethodDef CallMethods[] = {
  { "mgcv_tmm",(DL_FUNC)&mgcv_tmm,5}, 
  {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
    {"gsl_bspline",       (DL_FUNC) &gsl_bspline,        9},
    {"gsl_bspline_deriv", (DL_FUNC) &gsl_bspline_deriv, 11},
    {"RuniqueCombs",      (DL_FUNC) &RuniqueCombs,       4},
    {NULL, NULL, 0}
};

void R_init_ma(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
