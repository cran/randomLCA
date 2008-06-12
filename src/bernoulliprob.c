#include <R.h>
#include <Rinternals.h>

SEXP bernoulliprob(SEXP patterns, SEXP outcomep)
{
	SEXP ans;
	int irow, outcome, index, noutcomes, nrows;
	double *rpatterns = REAL(patterns), *routcomep = REAL(outcomep), *rans;
	double product;
	
	noutcomes = LENGTH(outcomep);
	nrows = LENGTH(patterns)/noutcomes;
	
	PROTECT(ans = allocVector(REALSXP,nrows));
	
	rans = REAL(ans);
	
	for (irow=0; irow < nrows; irow++) {
		product=1;
		for (outcome=0; outcome <noutcomes; outcome++) {
			index = irow+outcome*nrows;
			if (!ISNAN(rpatterns[index]))
				product = product*(rpatterns[index]*routcomep[outcome]+
					(1-rpatterns[index])*(1-routcomep[outcome]));
		}
		rans[irow]=product;
	}
	UNPROTECT(1);		
	return ans;
}
