#include <R.h>
#include <Rinternals.h>

SEXP lcemalgorithm(SEXP patterns, SEXP outcomep, SEXP classp, SEXP freq, SEXP verbose)
{
	SEXP ans, ill, ill2, classprob, logl, newoutcomep, newclassp;
	int emit, irow, iclass, outcome, index, noutcomes, nrows, nclass, lverbose;
	double *rpatterns = REAL(patterns),  
		*rfreq=REAL(freq), oldll, *rill, *rill2, *rclassprob, product,
		sumll, ll, sumfreq, sum1, sum2, *rlogl, *rnewoutcomep, *rnewclassp;

	lverbose = asLogical(verbose);

	nclass = LENGTH(classp);
	noutcomes = LENGTH(outcomep)/nclass;
	nrows = LENGTH(patterns)/noutcomes;
	
	PROTECT(ans = allocVector(VECSXP,3));

	PROTECT(logl = allocVector(REALSXP,1));
	PROTECT(ill = allocMatrix(REALSXP, nrows, nclass));
	PROTECT(ill2 = allocVector(REALSXP,nrows));
	PROTECT(classprob = allocMatrix(REALSXP,  nrows, nclass));
	PROTECT(newoutcomep = duplicate(outcomep));
	PROTECT(newclassp = duplicate(classp));
	
	rlogl = REAL(logl);
	rill = REAL(ill);
	rill2 = REAL(ill2);
	rclassprob = REAL(classprob);
	rnewoutcomep = REAL(newoutcomep);
	rnewclassp = REAL(newclassp);
	
	
	oldll = NA_REAL;
	emit =  0;
	while (TRUE) {
		/* calculate probabilities under each class */
		for (iclass=0; iclass < nclass; iclass++) {
			/* calculate the outcome probabilities for this class */
			for (irow=0; irow < nrows; irow++) {
				product=1;
				for (outcome=0; outcome < noutcomes; outcome++) {
					index = irow+outcome*nrows;
					if (!ISNAN(rpatterns[index]))
						product = product*(rpatterns[index]*rnewoutcomep[outcome*nclass+iclass]+
							(1.0-rpatterns[index])*(1.0-rnewoutcomep[outcome*nclass+iclass]));
				}
				rill[irow+iclass*nrows]=product*rnewclassp[iclass];
			}
		}
		ll =0.0;
		for (irow=0; irow < nrows; irow++) {
			sumll = 0.0;
			for (iclass=0; iclass < nclass; iclass++) {
				if (!ISNAN(rill[irow+iclass*nrows]))
					sumll = sumll+rill[irow+iclass*nrows];
			}
			rill2[irow] = sumll;
			if (!ISNAN(sumll))
				ll = ll+log(sumll)*rfreq[irow];
		}
	
		if (ISNAN(oldll)) oldll = 2*ll;
		if (fabs((ll-oldll)/ll) < 1.0e-8) break;
		oldll = ll;
		
		/* estimated posterior probability */
		for (irow=0; irow < nrows; irow++) {
			for (iclass=0; iclass < nclass; iclass++) {
				rclassprob[irow+iclass*nrows] = rill[irow+iclass*nrows]/rill2[irow];
			}
		}

		
		/* new rnewclassp */
		for (iclass=0; iclass < nclass; iclass++) rnewclassp[iclass]=0.0;
		sumfreq=0.0;
		for (irow=0; irow < nrows; irow++) {
			sumfreq += rfreq[irow];
			for (iclass=0; iclass < nclass; iclass++) 
				rnewclassp[iclass] += rclassprob[irow+iclass*nrows]*rfreq[irow];
		}
		for (iclass=0; iclass < nclass; iclass++) rnewclassp[iclass]/= sumfreq;

		
		/* new outcome probabilities */
		for (iclass=0; iclass < nclass; iclass++) {
			for (outcome=0; outcome < noutcomes; outcome++) {
				sum1=0;
				sum2=0;
				for (irow=0; irow < nrows; irow++) {
						if (!ISNAN(rpatterns[irow+nrows*outcome])) {
							sum1 += rpatterns[irow+nrows*outcome]*
								rclassprob[irow+nrows*iclass]*rfreq[irow];
							sum2 += rclassprob[irow+nrows*iclass]*rfreq[irow];
					}					
				}
				rnewoutcomep[iclass+nclass*outcome]= sum1/sum2;
			}
		}
		for (index=0; index < nclass*noutcomes; index++) {
			if (rnewoutcomep[index] < 1.0e-10)  rnewoutcomep[index] = 1.0e-10;
			if (rnewoutcomep[index] > 1-1.0e-10)  rnewoutcomep[index] = 1-1.0e-10;
		}
		emit++;
		if (((emit % 100)==0) & lverbose) Rprintf("iteration %d logl %f\n",emit,ll);
		if (emit > 250) break;
	}

	rlogl[0]=ll;	
	SET_VECTOR_ELT(ans, 0, logl);
	SET_VECTOR_ELT(ans, 1, newoutcomep);
	SET_VECTOR_ELT(ans, 2, newclassp);
	

	UNPROTECT(7);		
	return ans;
}
