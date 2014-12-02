#include <mex.h>
#include "tools.h"
#include "tv_core.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double delta, L, mu, eps;
	register double *y,*x,*epsilon_kf,*kf;
	mxArray *Ym;
	mxArray *zp;
	register int maxiter;  
	int type,m,n;

	if(nrhs != 7)
		printf("Should contain 7 input parameters but has %i\n",nrhs); DRAW

	Ym = (mxArray*) prhs[0]; /* Pointer to matrix structure：B*/ 
	y = mxGetPr(Ym); /* Pointer to the first element of matrix data B;(double type only)*/
	
	zp = (mxArray*) prhs[1];
	delta = (double)(mxGetScalar(zp));//pointer to first real element of matrix mxArray zp

	zp = (mxArray*) prhs[2];
	eps = (double)(mxGetScalar(zp));//epsilon

	zp = (mxArray*) prhs[3];
	L = (double)(mxGetScalar(zp));//Lmu

	zp = (mxArray*) prhs[4];
	mu = (double)(mxGetScalar(zp));//mu

	zp = (mxArray*) prhs[5];
	maxiter = (int)(mxGetScalar(zp)); // max iteration number

	zp = (mxArray*) prhs[6];
	type = (int)(mxGetScalar(zp)); // type: unknown

	m = mxGetM(Ym), n = mxGetN(Ym); // number of rows and column respectively

	/*Allocate memory and assign output pointer*/
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); /*mxReal is our data-type:X*/
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	/* Get a pointer to the data space in our newly allocated memory */
	x = mxGetPr(plhs[0]);//X
	kf = mxGetPr(plhs[1]);//k
	epsilon_kf = mxGetPr(plhs[2]);//epsilon_k

	if(type==1){
		tv_denoise_core(x,y,delta,eps,L,mu,m,n,maxiter,kf,epsilon_kf);
	}
	else{
		tv_denoise_core_org(x,y,delta,eps,L,mu,m,n,maxiter,kf,epsilon_kf);
	}

}
