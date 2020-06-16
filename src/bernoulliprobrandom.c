#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP bernoulliprobrandom(SEXP patterns, SEXP outcomex,SEXP lambdacoef, 
	SEXP gh, SEXP momentdata, SEXP probit)
{
	SEXP ans, temp;
	int irow, outcome, index, noutcomes, nrows, ipoint, npoints, level2size, ilambda, lprobit, *rpatterns = INTEGER(patterns);
	double  *routcomex = REAL(outcomex), *rans, *rtemp,
		neww,newp, *rmomentdata=REAL(momentdata),
		*rgh=REAL(gh),*rlambdacoef=REAL(lambdacoef);
	double onesum, sum, myoutcomex, myoutcomep, maxtemp;

	lprobit = asLogical(probit);
	
	noutcomes = LENGTH(outcomex);
	nrows = LENGTH(patterns)/noutcomes;
	npoints = LENGTH(gh)/2;
	level2size=LENGTH(lambdacoef);
	
	PROTECT(ans = allocVector(REALSXP,nrows));
	
	rans = REAL(ans);
	
	PROTECT(temp = allocVector(REALSXP,npoints));

	 rtemp = REAL(temp);	
	
	for (irow=0; irow < nrows; irow++) {
		/* Rprintf("irow  %d\n",irow); */
		sum=0.0;
/* calculate transformed w and p */
		for (ipoint=0; ipoint < npoints; ipoint++) {
			/* Rprintf("momentdata  %f,%f\n",rmomentdata[irow],rmomentdata[nrows+irow]); */
			newp = rmomentdata[irow]+rmomentdata[nrows+irow]*rgh[ipoint];
			neww = log(rmomentdata[nrows+irow])+
				(rgh[ipoint]*rgh[ipoint])/2.0+log(rgh[npoints+ipoint])-
				newp*newp/2.0;
			/* Rprintf("newp,neww  %f,%f\n",newp,neww); */
			ilambda=0;
			onesum=0.0;
			for (outcome=0; outcome <noutcomes; outcome++) {
				/* calculate outcome probability for this outcome */
				myoutcomex = routcomex[outcome]+
					rlambdacoef[ilambda]*newp;
				if (lprobit)
					myoutcomep=pnorm(myoutcomex,0,1,TRUE,TRUE);
				else
					myoutcomep=-log(1+exp(-myoutcomex));
				ilambda=(ilambda+1) % level2size;				
				/* update likelihood for this observation */
			/*  Rprintf("myoutcomep  %f\n",myoutcomep); */
				index = irow+outcome*nrows;
				if (rpatterns[index]!=NA_INTEGER) {
				  if (rpatterns[index]==1) onesum = onesum+myoutcomep;
				  else onesum = onesum+log(1-exp(myoutcomep)); 
				}
			}
			rtemp[ipoint] = onesum+neww;
		}
		maxtemp = R_NegInf;
		for (ipoint=0; ipoint < npoints; ipoint++)
		  if(rtemp[ipoint] > maxtemp) maxtemp = rtemp[ipoint];
		sum = 0.0;
		for (ipoint=0; ipoint < npoints; ipoint++)
		  sum = sum+exp(rtemp[ipoint]-maxtemp);
		rans[irow]=maxtemp+log(sum);
	}
	UNPROTECT(2);		
	return ans;
}
