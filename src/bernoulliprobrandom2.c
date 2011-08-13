#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP bernoulliprobrandom2(SEXP patterns, SEXP outcomex,SEXP lambdacoef, SEXP ltaucoef,
	SEXP gh, SEXP momentdata, SEXP probit, SEXP updatemoments)
{
	SEXP ans, ill, outcomep, newmomentdata, ll2,  ll3, sum3, sum4, sum5, sum3total, sum4total, sum5total, mye_2, mye2_2, mye23_2,e_2,e2_2,e23_2,e,e2;
	int irow, outcome, index, noutcomes, nrows, i2point, i3point, npoints, blocksize,
		ilambda, lprobit, lupdatemoments, nlevel2, i2;
	double *rpatterns = REAL(patterns), *routcomex = REAL(outcomex), *rans,
		*routcomep,new3w, new3p, new2w, new2p, *rmomentdata=REAL(momentdata),
		*rgh=REAL(gh),*rlambdacoef=REAL(lambdacoef), *rill, *rnewmomentdata,
		 *rll2,  *rll3, *rmye_2, *rmye2_2,  *rmye23_2, *re_2, *re2_2, *re23_2, *re, *re2,
		 product, sum2, sum2total, suml3, suml4,suml5, *rsum3, *rsum4, *rsum5, *rsum3total,
		 *rsum4total, *rsum5total, sumll3, llsum, myoutcomex, myoutcomep, rtaucoef;
	
	lprobit = asLogical(probit);
	lupdatemoments = asLogical(updatemoments);
	
	rtaucoef = exp(REAL(ltaucoef)[0]);
	
	noutcomes = LENGTH(outcomex);
	nrows = LENGTH(patterns)/noutcomes;
	npoints = LENGTH(gh)/2;
	blocksize=LENGTH(lambdacoef);
	nlevel2=noutcomes/blocksize;

	
	if (lupdatemoments) PROTECT(ans = allocVector(VECSXP,2));
	else PROTECT(ans = allocVector(VECSXP,1));

	PROTECT(ill = allocVector(REALSXP,nrows));
	rill = REAL(ill);

	if (lupdatemoments) {
		PROTECT(newmomentdata = allocMatrix(REALSXP,nrows,2+3*nlevel2));
		rnewmomentdata = REAL(newmomentdata);

		PROTECT(sum3 = allocVector(REALSXP,nlevel2));
		rsum3 = REAL(sum3);
		PROTECT(sum4 = allocVector(REALSXP,nlevel2));
		rsum4 = REAL(sum4);
		PROTECT(sum5 = allocVector(REALSXP,nlevel2));
		rsum5 = REAL(sum5);

		PROTECT(sum3total = allocVector(REALSXP,nlevel2));
		rsum3total = REAL(sum3total);
		PROTECT(sum4total = allocVector(REALSXP,nlevel2));
		rsum4total = REAL(sum4total);
		PROTECT(sum5total = allocVector(REALSXP,nlevel2));
		rsum5total = REAL(sum5total);		
	}


	for (irow=0; irow < nrows; irow++) {
		/* Rprintf("irow  %d\n",irow); */
		suml3=0.0;
		if (lupdatemoments) {
			suml4=0.0;
			suml5=0.0;
			for (i2=0; i2<nlevel2; i2++) {			
				rsum3total[i2]=0;
				rsum4total[i2]=0;
				rsum5total[i2]=0;
			}
		}
		for (i3point=0; i3point < npoints; i3point++) {
			/* Rprintf("momentdata  %f,%f\n",rmomentdata[irow],rmomentdata[nrows+irow]); */
			new3p = rmomentdata[irow]+rmomentdata[nrows+irow]*rgh[i3point];
			new3w = log(sqrt(2.0*M_PI))+log(rmomentdata[nrows+irow])+
				(rgh[i3point]*rgh[i3point])/2.0+log(rgh[npoints+i3point])+
				dnorm(new3p,0.0,1.0,TRUE);
			sum2total=0;
			for (i2=0; i2<nlevel2; i2++) {			
				sum2=0.0;
				if (lupdatemoments) {
				/* needs to be repeated for each level 2 unit */
					rsum3[i2]=0.0;
					rsum4[i2]=0.0;
					rsum5[i2]=0.0;
				}
				for (i2point=0; i2point < npoints; i2point++) {
					/* calculate quadrature points for level2 unit */
					/* Rprintf("momentdata  %f,%f\n",rmomentdata[irow+2+3*i2],rmomentdata[nrows+irow+3+3*i2]); */
					new2p = rmomentdata[irow+nrows*(2+3*i2)]+rmomentdata[irow+nrows*(3+3*i2)]*rgh[i2point];
					new2p = new2p-rmomentdata[irow+nrows*(4+3*i2)]*(new3p-rmomentdata[irow]);
					new2w = log(sqrt(2.0*M_PI))+log(rmomentdata[irow+nrows*(3+3*i2)])+
						(rgh[i2point]*rgh[i2point])/2.0+log(rgh[npoints+i2point])+
						dnorm(new2p,0.0,1.0,TRUE);
					/* calculate logl for level 2 unit */
					product=1.0;
					for (outcome=0; outcome <blocksize; outcome++) {
						/* calculate outcome probability for this outcome */
						myoutcomex = routcomex[outcome+i2*blocksize]+ (new3p+new2p*rtaucoef)*rlambdacoef[outcome];
						if (lprobit)
							myoutcomep=pnorm(myoutcomex,0,1,TRUE,FALSE);
						else
							myoutcomep=exp(myoutcomex)/(1+exp(myoutcomex));
						/* update likelihood for this observation */
					/*  Rprintf("myoutcomep  %f\n",myoutcomep); */
						index = irow+(outcome+i2*blocksize)*nrows;
						if (!ISNAN(rpatterns[index])) {
							product = product*(rpatterns[index]*myoutcomep+
								(1-rpatterns[index])*(1-myoutcomep));
						}
					}
					/* sum across level 2 quadrature points */
					sum2=sum2+product*exp(new2w);
					if (lupdatemoments) {
						rsum3[i2]=rsum3[i2]+new2p*product*exp(new2w);
						rsum4[i2]=rsum4[i2]+new2p*new2p*product*exp(new2w);
						rsum5[i2]=rsum5[i2]+new3p*new2p*product*exp(new2w);
						/*
						if (irow==0)  {
							Rprintf("i2point %d\n",i2point);
							Rprintf("sum3  %f\n",rsum3[i2]);
							Rprintf("sum4  %f\n",rsum4[i2]);
							Rprintf("sum5  %f\n",rsum5[i2]);
						}
						*/
					}
				}
				/* sum across level2 units */
				sum2total=sum2total+log(sum2);
				if (lupdatemoments) {
					rsum3[i2]=rsum3[i2]/sum2;
					rsum4[i2]=rsum4[i2]/sum2;
					rsum5[i2]=rsum5[i2]/sum2;
					/*
					if (irow==0)  {
						Rprintf("sum2 %f\n",sum2);
						Rprintf("sum3  %f\n",rsum3[i2]);
						Rprintf("sum4  %f\n",rsum4[i2]);
						Rprintf("sum5  %f\n",rsum5[i2]);
					}
					*/
				}
			}
			/* if (irow==0)  Rprintf("level2 ll  %f\n",sum2total); */
			suml3=suml3+exp(sum2total+new3w);
			if (lupdatemoments) {
				suml4=suml4+new3p*exp(sum2total+new3w);
				suml5=suml5+new3p*new3p*exp(sum2total+new3w);
				for (i2=0; i2<nlevel2 ;i2++) {
					rsum3total[i2]=rsum3total[i2]+exp(sum2total+new3w)*rsum3[i2];
					rsum4total[i2]=rsum4total[i2]+exp(sum2total+new3w)*rsum4[i2];
					rsum5total[i2]=rsum5total[i2]+exp(sum2total+new3w)*rsum5[i2];
					/*
						if (irow==0)  {
							Rprintf("sum3total  %f\n",rsum3total[i2]);
							Rprintf("sum4total  %f\n",rsum4total[i2]);
							Rprintf("sum5total  %f\n",rsum5total[i2]);
						}
					*/
				}
			}
			/* if (irow==0)  Rprintf("new3w  %f\n",new3w); */
			/* if (irow==0)  Rprintf("suml3  %f\n",suml3); */
		}
		rill[irow]=log(suml3);
		if (lupdatemoments) {
			rnewmomentdata[irow]=suml4/suml3;
			rnewmomentdata[nrows+irow]=sqrt((suml5-suml4*suml4/suml3)/suml3);
			/*		if (irow==0)  {
						Rprintf("e  %f\n",rnewmomentdata[irow]);
						Rprintf("e2  %f\n",rnewmomentdata[nrows+irow]);
					}
			*/
		}
		/* if (irow==0)  Rprintf("level3 ll  %f\n",rill[irow]); */
		if (lupdatemoments) {
			for (i2=0; i2<nlevel2; i2++) {
				rnewmomentdata[(3*i2+2)*nrows+irow]=rsum3total[i2]/suml3;
				rnewmomentdata[(3*i2+3)*nrows+irow]=sqrt(rsum4total[i2]/suml3-
					rnewmomentdata[(3*i2+2)*nrows+irow]*rnewmomentdata[(3*i2+2)*nrows+irow]);
				rnewmomentdata[(3*i2+4)*nrows+irow]=-(rsum5total[i2]/suml3-
					rnewmomentdata[irow]*rnewmomentdata[(3*i2+2)*nrows+irow])/
					(rnewmomentdata[nrows+irow]*rnewmomentdata[nrows+irow]);
					/*
					if (irow==0)  {
						Rprintf("e_2  %f\n",rnewmomentdata[(3*i2+2)*nrows+irow]);
						Rprintf("e_22  %f\n",rnewmomentdata[(3*i2+3)*nrows+irow]);
						Rprintf("e_23  %f\n",rnewmomentdata[(3*i2+4)*nrows+irow]);
					}
					*/
			}
		}
	}
	SET_VECTOR_ELT(ans, 0, ill);
	if (lupdatemoments) SET_VECTOR_ELT(ans, 1, newmomentdata);
	/* how many */
	if (lupdatemoments) UNPROTECT(9);
	else UNPROTECT(2);
	return ans;
}
