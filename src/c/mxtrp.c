#include <mex.h>
#include "tools.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double delta, nu1v=0;
	register double *x,*p,*c,*t,*nu1,*cs;
	mxArray *cm,*pm;
	mxArray *zp;
	int m,n,mn;  

	/*
    Solves the TRSP
    
    minimize    x'*P*x - 2*c'*x
    subject to  x'*x <= delta

    where P := diag(p).
	*/

	if(nrhs != 3)
		printf("Should contain 3 input parameters but has %i\n",nrhs);

	cm = (mxArray*)prhs[0]; /* Pointer to matrix structure*/
	c = mxGetPr(cm); /* Pointer to the matrix data*/

	pm = (mxArray*)prhs[1];
	p = mxGetPr(pm);

	zp = (mxArray*)prhs[2];
	delta = (double)(mxGetScalar(zp));

  m = mxGetM(cm); n = mxGetN(cm);

	/*Allocate memory and assign output pointer*/
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); /*mxReal is our data-type*/

	/* Get a pointer to the data space in our newly allocated memory */
	x = mxGetPr(plhs[0]);

	mn=m*n;
	cs = malloc(mn*sizeof(double));
	t  = malloc(mn*sizeof(double)); /* If trp2 is called multiple times, we can avoid multiple mallocs for this temp vector */

  nu1 = &nu1v;

	trp2(c,p,delta,nu1,t,cs,mn,x); /* solution ends in x */

	free(t);free(cs);
} 
