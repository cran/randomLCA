#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP bernoulliprob(SEXP patterns, SEXP outcomep)
{
	SEXP ans;
	int irow, outcome, index, noutcomes, nrows, *rpatterns = INTEGER(patterns);
	double  *routcomep = REAL(outcomep), *rans;
	double product;
	
	noutcomes = LENGTH(outcomep);
	nrows = LENGTH(patterns)/noutcomes;
	
	PROTECT(ans = allocVector(REALSXP,nrows));
	
	rans = REAL(ans);
	
	for (irow=0; irow < nrows; irow++) {
		product=1;
		for (outcome=0; outcome <noutcomes; outcome++) {
			index = irow+outcome*nrows;
		  if (rpatterns[index]!=NA_INTEGER) {
		    if (rpatterns[index]==1) product = product*routcomep[outcome];
		    else product = product*(1-routcomep[outcome]); 
		  }
		}
		rans[irow]=product;
	}
	UNPROTECT(1);		
	return ans;
}
