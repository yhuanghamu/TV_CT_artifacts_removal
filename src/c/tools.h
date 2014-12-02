#ifndef __TOOLS__
#define __TOOLS__

#if (defined(_WIN32) || defined(__WIN32__) )
#include <mex.h>
#define DRAW mexEvalString("drawnow;");
#else
#define DRAW ;
#endif

#ifdef FFTW3
#include <fftw3.h>
#endif

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))

void trp(double *c,double *d,double *dis,double nuc,double *y,double delta,double *nu_old,double *t,double *qs,double *q,int mn,double maxy,double mind,double *x);

void trp2(double *c,double *p,double delta,double *nu_old,double *t,double *cs,int mn,double *x);

double minf(double *x,int N);

#ifdef FFTW3
void odct2(double *x,int m,int n, fftw_plan *p,int k);

void oidct2(double *x,int m, int n, fftw_plan *p, int k);
#endif

double ddot_(int *n, double *x, int *incx, double *y, int *incy);

double dnrm2_(int *n, double *x, int *incx);

void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
#endif

//daxpy:compute y := alpha * x + y  DAXPY(N, ALPHA, X, INCX, Y, INCY)
//	https://docs.oracle.com/cd/E19422-01/819-3691/daxpy.html
//dnrm2_: dnrm2_ := sqrt( x'*x )
//	http://www.netlib.org/lapack/explore-html/da/d7f/dnrm2_8f.html
//ddot_: ddot_ := x.*y forms the dot product of two vectors.
//	http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html

