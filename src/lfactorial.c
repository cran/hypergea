
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>


void intern_lfactorial( int *xx, int *L, double *r ){
	int l=0;
	int LL=(*L);
	double *preCalcFact;
	int mx=xx[0];
	for( l=1; l<LL; l++ ){
		if((xx[l])>mx){ mx = xx[l]; }
	}
	if( mx%2==1 ){mx++;}

	preCalcFact=(double*)malloc((mx+1)*sizeof(double));
	preCalcFact[0]=0;
	for( l=1; l<= mx; l++){
		preCalcFact[l]=preCalcFact[l-1]+log((float)l);
		l++;
		preCalcFact[l]=preCalcFact[l-1]+log((float)l);
	}

	*r=preCalcFact[(int)xx[0]];
	for(l=1; l<LL; l++){
		*r += preCalcFact[(int)xx[l]];
	}

	free(preCalcFact);
}
SEXP lfactorial(SEXP x){
	int *xx;
	int LL=0;
	double r;
	SEXP LM;

	xx=INTEGER(x); LL=length(x);
	intern_lfactorial(xx, &LL, &r);
	LM=PROTECT(allocVector(REALSXP, 1));
	REAL(LM)[0]=r;
	UNPROTECT(1);
return(LM);
}

