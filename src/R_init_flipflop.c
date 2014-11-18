#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "Rinterface.h"

#include "R_init_flipflop.h"

R_CMethodDef cMethods[] = {
  //  {"ffProcesssam", (DL_FUNC) & ffProcesssam, 6, {STRSXP, STRSXP, STRSXP, STRSXP, STRSXP, STRSXP}},
  {"ffProcesssam", (DL_FUNC) & ffProcesssam, 8},
  {NULL, NULL, 0}
};

R_CallMethodDef callMethods[] = {
  /*
  {"R_SWIG_debug_getCallbackFunctionData", 
   (DL_FUNC) & R_SWIG_debug_getCallbackFunctionData, 0},
  {"R_SWIG_R_pushCallbackFunctionData", 
   (DL_FUNC) & R_SWIG_R_pushCallbackFunctionData, 2},
  */
  {"R_swig_multLeftDiag",
   (DL_FUNC) & R_swig_multLeftDiag, 3},
  {"R_swig_fistaFlat",
   (DL_FUNC) & R_swig_fistaFlat, 39},
  {"R_swig_evalPathCoding",
   (DL_FUNC) & R_swig_evalPathCoding, 21},
  {"R_swig_sepCostsPathCoding",
   (DL_FUNC) & R_swig_sepCostsPathCoding, 19},
  {"R_swig_solverPoisson",
   (DL_FUNC) & R_swig_solverPoisson, 8},
  {"R_swig_solverPoissonFull",
   (DL_FUNC) & R_swig_solverPoissonFull, 8},
  {NULL, NULL, 0}
};


void R_init_flipflop(DllInfo * info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

void R_unload_flipflop(DllInfo *info)
{
    (void) info;
}
