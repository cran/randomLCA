#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

SEXP bernoulliprob(SEXP patterns, SEXP outcomep);

void fixedloglik(int n, double *pars, double *g, void *state)
{
	SEXP patterns, freq, mynrows, mynclass, mynoutcomes, probit, outcomep, onell,
		mytempstorage;
	int wantprobit, nrows, nclass, noutcomes, class, outcome, irow;
	double *classx, *classp,*outcomex, *routcomep, *ill2, *ronell, *rfreq,
		sum1,ll, *tempstorage;
/*	Rprintf("enetring loglik\n"); */
/* find each of the components of state */

	patterns = VECTOR_ELT((SEXP)state,0);
	freq = VECTOR_ELT((SEXP)state,1);
	rfreq=REAL(freq);
	mynclass = VECTOR_ELT((SEXP)state,2);
	nclass=asInteger(mynclass);	
	mynrows = VECTOR_ELT((SEXP)state,3);
	nrows=asInteger(mynrows);	
	mynoutcomes = VECTOR_ELT((SEXP)state,4);
	noutcomes=asInteger(mynoutcomes);	
	probit = VECTOR_ELT((SEXP)state,5);
	wantprobit=asLogical(probit);

	mytempstorage = VECTOR_ELT((SEXP)state,6);
	tempstorage=REAL(mytempstorage);
	
	 outcomep = VECTOR_ELT((SEXP)state,7);
    routcomep=REAL(outcomep);
/*
	PrintValue(patterns);
	PrintValue(freq);
	PrintValue(mynclass);
	PrintValue(mynrows);
	PrintValue(mynoutcomes);
	PrintValue(probit);
	PrintValue(mytempstorage);
	PrintValue(outcomep);

*/

/* extract the parameters */
/* class probabilities */
	classx = &tempstorage[0];
	classx[0]=0.0;
	for (class=1; class <nclass; class++) {
		classx[class]=pars[class-1];
	}
	classp = &tempstorage[nclass];
	sum1=0.0;
	for (class=0; class <nclass; class++) {
		sum1=sum1+exp(classx[class]);
	}
	for (class=0; class <nclass; class++) {
		classp[class]=exp(classx[class])/sum1;
	}
/*	for (class=0; class <nclass; class++) {
		Rprintf(" %f",classp[class]);
	}
	Rprintf("\n");
*/
/* outcome probabilities */
	outcomex = &tempstorage[2*nclass];
	for (outcome=0; outcome <noutcomes*nclass; outcome++) {
		outcomex[outcome]=pars[nclass+outcome-1];
	}
/* storage for ill */
	ill2 = &tempstorage[2*nclass+noutcomes*nclass];
	for (class=0; class <nclass; class++) {
		for (outcome=0; outcome <noutcomes; outcome++) {
			if (wantprobit) {
				routcomep[outcome] = pnorm(outcomex[class+outcome*nclass],0,1,TRUE,FALSE);
			}
			else {
				routcomep[outcome] = 1/(1+exp(-outcomex[class+outcome*nclass]));
			}
/*			Rprintf(" %f",routcomep[outcome]); */
		}
/*		Rprintf("\n");
	PrintValue(outcomep);
*/
	PROTECT(onell = bernoulliprob(patterns,outcomep));
		ronell=REAL(onell);
		if (class==0) {
			for (irow=0;irow <nrows;irow++) {
				ill2[irow]=ronell[irow]*classp[class];
			}
		} else {
			for (irow=0;irow <nrows;irow++) {
				ill2[irow]=ill2[irow]+ronell[irow]*classp[class];
			} 
		}
		UNPROTECT(1);
	}
/* calculate total loglikelihood */
	ll = 0.0;
	for (irow=0;irow <nrows;irow++) {
		ll = ll -log(ill2[irow])*rfreq[irow];
	} 
/*	Rprintf(" ll = %f\n",ll); */
	*g=ll;
/*	Rprintf("exitting loglik\n"); */
}

static void d1fcn_dum(int n, double *x, double *g, void *state)
{
/*	dummy routine to prevent unsatisfied external diagnostic
 *	when specific analytic gradient function not supplied. */
}

static void d2fcn_dum(int nr, int n, double *x, double *h, void *state)
{
/*	dummy routine to prevent unsatisfied external diagnostic
 *	when specific analytic hessian function not supplied. */
}

SEXP fixednlm(SEXP p, SEXP state, SEXP dohessian, SEXP printlevel, SEXP gradtl, SEXP itnlim)
{

    SEXP ans, xpls, gpls,  a,  myfpls, mycode, myiterations;

    double *rp, *typsiz, fscale, *rgradtl, stepmx, steptl,
	steptol, *rxpls, *rgpls,  *ra, *wrk, dlt, ll, *rlogl, fpls, *xfpls;

    int icode, i, j, k,  method, iexp, msg,
	nparams, ndigit,  want_hessian, itncnt, 
	*rprintlevel, *ritnlim, *xcode, *xiterations;

	/* get parameters */
	
	rp = REAL(p);
	nparams = LENGTH(p);
	
	want_hessian = asLogical(dohessian);

	rprintlevel = INTEGER(printlevel);
	rgradtl = REAL(gradtl);
	ritnlim = INTEGER(itnlim);

	/* set up options */
	
	typsiz = (double*)R_alloc(nparams, sizeof(double));
	for (i = 0; i < nparams; i++)
		typsiz[i] = 1.0;
	fscale = 1.0;
    method = 1;	/* Line Search */
    iexp =  1; /* Function calls are not expensive */
	if (*rprintlevel==2) msg = 25;
	else msg=9;
	
	ndigit = 12;
    dlt = 1.0;
    /* rgradtl from inputs */
    stepmx = 1000.0;
    steptl=1.0e-6;

	PROTECT(ans = allocVector(VECSXP,6));
    
    PROTECT(xpls = allocVector(REALSXP,nparams));
    rxpls=REAL(xpls);
    
    PROTECT(gpls = allocVector(REALSXP,nparams));
    rgpls=REAL(gpls);

    PROTECT(a = allocMatrix(REALSXP,nparams,nparams));
    ra=REAL(a);

    wrk = (double*)R_alloc(8*nparams, sizeof(double));

/*	Rprintf("calling optif9\n");*/
	
	optif9(nparams, nparams, rp, (fcn_p) fixedloglik, (fcn_p)d1fcn_dum, (d2fcn_p)d2fcn_dum,
       state , typsiz, fscale, method,
       iexp, &msg, ndigit, *ritnlim, 0, 0,
       dlt, *rgradtl, stepmx,steptl,
       rxpls, &fpls, rgpls, &icode, ra,
       wrk, &itncnt);
/*	Rprintf("after optif9\n"); */

    if (msg < 0)
		error("warning optif9 (msg = %d)", msg);

/*
	ignore warnings these will be looked at in calling routine
	if (code != 0 && (omsg&8) == 0)
	optcode(code);
*/

    if (want_hessian) {
		fdhess(nparams, rxpls, fpls, (fcn_p) fixedloglik, state, ra, nparams, &wrk[0], &wrk[nparams],
	       ndigit, typsiz);
		for (i = 0; i < nparams; i++)
	    	for (j = 0; j < i; j++)
				ra[i + j * nparams] = ra[j + i * nparams];
	}
	SET_VECTOR_ELT(ans, 0, xpls);

	PROTECT(myfpls = allocVector(REALSXP,1));
	xfpls = NUMERIC_POINTER(myfpls);
	xfpls[0]=fpls;

	SET_VECTOR_ELT(ans, 1, myfpls);
	SET_VECTOR_ELT(ans, 2, gpls);
	if (want_hessian) SET_VECTOR_ELT(ans, 3, a);

	PROTECT(mycode = allocVector(INTSXP,1));
	xcode = INTEGER_POINTER(mycode);
	xcode[0]=icode;
	SET_VECTOR_ELT(ans, 4, mycode);

	PROTECT(myiterations = allocVector(INTSXP,1));
	xiterations = INTEGER_POINTER(myiterations);
	xiterations[0]=itncnt;
	SET_VECTOR_ELT(ans, 5, myiterations);


/*
	PrintValue(xpls);
	PrintValue(myfpls);
	PrintValue(gpls);
	if (want_hessian) PrintValue(a);
	PrintValue(mycode);
*/
	UNPROTECT(7);
	
	return ans;
}
