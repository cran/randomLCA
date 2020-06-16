#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP bernoulliprob(SEXP patterns, SEXP loutcomep, SEXP lnoutcomep)
{
	SEXP ans;
	int irow, outcome, index, noutcomes, nrows, *rpatterns = INTEGER(patterns);
	double  *rloutcomep = REAL(loutcomep), *rlnoutcomep = REAL(lnoutcomep), *rans;
	double lproduct;
	
	noutcomes = LENGTH(loutcomep);
	nrows = LENGTH(patterns)/noutcomes;
	
	PROTECT(ans = allocVector(REALSXP,nrows));
	
	rans = REAL(ans);
	
	for (irow=0; irow < nrows; irow++) {
		lproduct=0;
		for (outcome=0; outcome <noutcomes; outcome++) {
			index = irow+outcome*nrows;
		  if (rpatterns[index]!=NA_INTEGER) {
		    if (rpatterns[index]==1) lproduct = lproduct+rloutcomep[outcome];
		    else lproduct = lproduct+rlnoutcomep[outcome]; 
		  }
		}
		rans[irow]=lproduct;
	}
	UNPROTECT(1);		
	return ans;
}
