#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include "tools.h"

/* Used at trp as the relatice accuracy needed*/
#define NTTOL_REL 1e-4

#ifdef BLAS

extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern double dnrm2_(int *n, double *x, int *incx);
extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);

#else

/* Only implements the simplest version of ddot which is used in this code 
   here, with incrx and incry = 1 hardcoded */
//  dobj = ddot_(&mn, df, &one, y, &one) - delta*dnrm2_(&mn, df, &one);
double ddot_(int *p_n, double *p_x, int *incx, double *p_y, int *incy)
{
     register int i, n;
     register double *x, *y;
     register double ddot = 0;
     x = p_x;
     y = p_y;
     n = *p_n;
     for (i = 0; i < n; ++i)
          ddot += (*x++)*(*y++);

     return ddot;
}

double dnrm2_(int *p_n, double *p_x, int *incx)
{
	register int i,n;
	register double *x;
	register double nrm2 = 0;
	
	x = p_x;
	n = *p_n;

	for (i = 0; i < n; ++i){
			nrm2 += (*x)*(*x);
			x++; 
	}

	return sqrt(nrm2);
}

/* Only implements the simplest version of daxpy which is used in this code 
   here, with incrx and incry = 1 hardcoded. *y points at the results (overwrite) */
void daxpy_(int *p_n, double *p_alpha, double *p_x, int *incx, double *p_y, int *incy)
{
	   register int i, n;
     register double *x, *y,alpha;
     x = p_x;
     y = p_y;
     n = *p_n;
		 alpha = *p_alpha;

     for (i = 0; i < n; ++i){
          *y = alpha*(*x) + (*y);
					y++;x++;
		 }
}

#endif

double minf(double *x,int N){
	double c;
	int i;

	c = x[0];

	for(i=1;i<N; i++)
		if(x[i]<c)
			c = x[i];
		
	return c;
}

#ifdef FFTW3

void odct2(double *x,int m,int n,	fftw_plan *p,int k){
  fftw_r2r_kind kind;
	int mn = m*n;
	int i;
	double c;

  kind = FFTW_REDFT10; /* type 2 dct2 */

	if(k == 0){
		*p = fftw_plan_r2r_2d(n, m, x, x, kind, kind, FFTW_ESTIMATE);
	}

  fftw_execute(*p);

	/* orthonalize the dct2 */
	c =sqrt(4.0/(mn))*0.25;

	for(i=0;i<mn;i++)
		x[i] *= c;

	c = sqrt(1.0/2.0);

	for(i=0;i<m;i++)
		x[i] *= c;

	for(i=0;i<n;i++)
		x[i*m] *= c;
}

void oidct2(double *x,int m, int n, fftw_plan *p,int k){

  fftw_r2r_kind kind;
	int mn = m*n;
	int i;
	double c;

  kind = FFTW_REDFT01; /* type 2 dct2 */

	/* orthonalize the dct2 */
	c = 1/(sqrt(m*n*4.0));   /* need other normalization because the idct is a scaled dct */
  		
	for(i=0;i<mn;i++)
		x[i] *= c;

	c = 1/sqrt(1.0/2.0);

	for(i=0;i<m;i++)
		x[i] *= c;

	for(i=0;i<n;i++)
    x[i*m] *= c;

	if(k == 0){
		*p = fftw_plan_r2r_2d(n, m, x, x, kind, kind, FFTW_ESTIMATE);
	}

  fftw_execute(*p);   
}
#endif
	
void trp(double *c,double *d,double *dis,double nuc,double *y,double delta,double *nu_old,double *t,double *qs,double *q,int mn,double nrm2y,double mind,double *x){
	/*
		Solves the trp problem
		
		min 1/2 x'x-c'x
		s.t. ||diag(d) x - y|| <= delta

		dis = inv(diag(d))**2,nuc=max(-dis),t=pointer to a temp vector,mn=m*n

    dis = d**-2
    nuc = max(-dis) 
	 */

	int one=1;
	register double tv,tvv,f,df,nu;
	register double NTTOL;

	register int i,k;

	
  NTTOL = nrm2y*NTTOL_REL*mind;

	NTTOL = MAX(1e-12,NTTOL);


	/* If c in constraint then this is the solution 
    if ||diag(d)c-y||_2<=delta */
	for (i=0; i<mn; i++){
		t[i] = d[i]*c[i]-y[i];
	
		/* Also calculate q and qs now we are running the loop*/
		q[i] = c[i]/d[i]-y[i]*dis[i];
		qs[i] = q[i]*q[i];
	}


	if(dnrm2_(&mn,t,&one) < delta){
		for (i=0; i<mn; i++){
			x[i] = c[i];
		}
		nu_old[0] = 0.0;
		return;
	}

	nu = nuc +1;

	if(nu<nu_old[0]) /* i.e. max function*/
		nu = nu_old[0];

	for(k=0;k < 30;k++){

		    /*f  =    sum( div(qs, (dis + nu)**2) ) - delta**2
					df = -2*sum( div(qs, (dis + nu)**3) ) */
		f=0;df=0;
		for (i=0; i<mn; i++){
			/*tv = dis[i]+nu;*/
			tv = dis[i]+nu;
			tvv = tv*tv;
			f += qs[i]/tvv;
			df += qs[i]/(tvv*tv);
		}
		f = f - delta*delta;
		df = -2*df;
			
    /*  if abs(f) < NTTOL: break */
		if(fabs(f)<NTTOL)
			break;

	/*	printf("f df %e %e\n",f,df); */

		/*	dnu = -f/df;
				a = 1.0;
				while(nu+a*dnu < nuc)
				   a *= 0.5;
								
			  nu += a*dnu; */
		tv = -f/df;
		tvv = 1;
		while(nu+tvv*tv < nuc)	
			tvv *= 0.5;
		
		nu += tvv*tv;
	}

	 if(k == 30){
		printf("Max number of iteration in trp reached. \n"); DRAW
		/* raise(SIGTERM) ; */
		}

	/*    x = div(y + div(q, dis + nu), d)
				e = blas.nrm2(mul(d,x)-y) */

	for (i=0; i<mn; i++){
		x[i] = (y[i]+( q[i]/(dis[i]+nu) ) ) / d[i];
	  t[i] =d[i]*x[i]-y[i];
	}

	if (dnrm2_(&mn,t,&one) > delta + sqrt(NTTOL)){
		printf("Inaccurate TRP solution. Proceed with care. \n"); DRAW
		return;
	}

	nu_old[0] = nu;
		
}

void trp2(double *c,double *p,double delta,double *nu_old,double *t,double *cs,int mn,double *x){
	/*
    Solves the TRSP
    
    minimize    x'*P*x - 2*c'*x
    subject to  x'*x <= delta

    where P := diag(p).
	*/

	int one=1;
	register double tv,tvv,f,df,nu,mp,nuc;
	register int i,k;


	mp = p[0];
	for(i=1;i<mn;i++)
		if(p[i]<mp)
			mp=p[i];
	
	if(mp>0.0){
		for(i=0;i<mn;i++)
			x[i] = c[i]/p[i];

		if(dnrm2_(&mn,x,&one) < delta){
			nu_old[0] = 0.0;
			return;			 /*return x */
			}
	}

	for(i=0;i<mn;i++)
			cs[i] = c[i]*c[i];

	nuc = -mp; /*max(-p)=-min(p)*/
							 
	nu = nuc +1;
							 
	/*if( nu < nu_old[0] ) i.e. max function
		nu = nu_old[0]; */


	for(k=0;k < 30;k++){

    /*    f  =    sum( div(cs, (p + nu)**2) ) - delta
          df = -2*sum( div(cs, (p + nu)**3) ) */

		f=0;df=0;
		for (i=0; i<mn; i++){
			tv = p[i]+nu;
			tvv = tv*tv;
			f += cs[i]/tvv;
			df += cs[i]/(tvv*tv);
		}
		f = f - delta;
		df = -2*df;
			
    /*  if abs(f) < NTTOL: break */
		/* if(abs(f)<NTTOL)
			 break; */

		/*	dnu = -f/df;
				a = 1.0;
				while(nu+a*dnu < nuc)
				   a *= 0.5;
								
			  nu += a*dnu; */
		tv = -f/df;
		tvv = 1;
		while(nu+tvv*tv < nuc)	
			tvv *= 0.5;
		
		nu += tvv*tv;
	}

	/*if(k == 100){
		printf("Max number of iteration in trp reached. \n");
		raise(SIGTERM) ;
		}*/

  /*  x = div(c, p + nu)
      e = delta - blas.nrm2(x)**2 */

	for (i=0; i<mn; i++){
		x[i] = c[i]/(p[i] + nu);
	}
	
	tv = dnrm2_(&mn,x,&one);
	tvv = tv*tv;

	if ( delta - tvv > 1e-7){
		printf("Inaccurate TRP solution. Proceed with care. \n"); DRAW
		return;
	}

	nu_old[0] = nu;
		
}
