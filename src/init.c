#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP lcemalgorithm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bernoulliprob(SEXP, SEXP, SEXP);
extern SEXP bernoulliprobrandom(SEXP, SEXP, SEXP,  SEXP, SEXP, SEXP);
extern SEXP bernoulliprobrandom2(SEXP, SEXP, SEXP, SEXP,
	SEXP, SEXP, SEXP, SEXP, SEXP);
	
static const R_CallMethodDef CallEntries[] = {
    {"lcemalgorithm",                        (DL_FUNC) &lcemalgorithm,                         7},
    {"bernoulliprob",                        (DL_FUNC) &bernoulliprob,                         3},
    {"bernoulliprobrandom",                  (DL_FUNC) &bernoulliprobrandom,                   6},
    {"bernoulliprobrandom2",                 (DL_FUNC) &bernoulliprobrandom2,                  9},
    {NULL, NULL, 0}
};

void R_init_randomLCA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
