#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP bernoulliprobrandom(SEXP patterns, SEXP outcomex,SEXP lambdacoef, SEXP gh,
	SEXP momentdata, SEXP probit)
{
	SEXP ans,outcomep;
	int irow, outcome, index, noutcomes, nrows, ipoint, npoints, blocksize, ilambda, isprobit;
	double *rpatterns = REAL(patterns), *routcomex = REAL(outcomex), *rans,
		*routcomep,neww,newp, *rmomentdata=REAL(momentdata),
		*rgh=REAL(gh),*rlambdacoef=REAL(lambdacoef);
	double product, sum, myoutcomex, myoutcomep;
	
	isprobit = asLogical(probit);
	
	noutcomes = LENGTH(outcomex);
	nrows = LENGTH(patterns)/noutcomes;
	npoints = LENGTH(gh)/2;
	blocksize=LENGTH(lambdacoef);
	
	PROTECT(ans = allocVector(REALSXP,nrows));
	
	rans = REAL(ans);
	
	
	for (irow=0; irow < nrows; irow++) {
		/* Rprintf("irow  %d\n",irow); */
		sum=0.0;
/* calculate transformed w and p */
		for (ipoint=0; ipoint < npoints; ipoint++) {
			/* Rprintf("momentdata  %f,%f\n",rmomentdata[irow],rmomentdata[nrows+irow]); */
			newp = rmomentdata[irow]+rmomentdata[nrows+irow]*rgh[ipoint];
			neww = log(sqrt(2.0*M_PI))+log(rmomentdata[nrows+irow])+
				(rgh[ipoint]*rgh[ipoint])/2.0+log(rgh[npoints+ipoint])+
				dnorm(newp,0.0,1.0,TRUE);
			/* Rprintf("newp,neww  %f,%f\n",newp,neww); */
			ilambda=0;
			product=1.0;
			for (outcome=0; outcome <noutcomes; outcome++) {
				/* calculate outcome probability for this outcome */
				myoutcomex = routcomex[outcome]+
					rlambdacoef[ilambda]*newp;
				if (isprobit)
					myoutcomep=pnorm(myoutcomex,0,1,TRUE,FALSE);
				else
					myoutcomep=exp(myoutcomex)/(1+exp(myoutcomex));
				ilambda=(ilambda+1) % blocksize;				
				/* update likelihood for this observation */
			/*  Rprintf("myoutcomep  %f\n",myoutcomep); */
				index = irow+outcome*nrows;
				if (!ISNAN(rpatterns[index])) {
					product = product*(rpatterns[index]*myoutcomep+
						(1-rpatterns[index])*(1-myoutcomep));
				}
			}
			sum=sum+product*exp(neww);
		}
		rans[irow]=sum;
	}
	UNPROTECT(1);		
	return ans;
}
