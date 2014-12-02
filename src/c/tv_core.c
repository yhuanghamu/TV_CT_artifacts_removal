#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tools.h"

/* Settings which makes the user do a CTRL-C break out of the loop*/
#if defined(LIBUT) && (defined(_WIN32) || defined(__WIN32__) )

#define STOPMARK utIsInterruptPending()
#define INITBREAK ;
bool utIsInterruptPending(void);

#else

#include <signal.h>
#define INITBREAK   sigint_cont = 1;  (void) signal(SIGINT , ex_sigint);
#define STOPMARK sigint_cont==0
int sigint_cont = 1;
void ex_sigint(int sig) {
	sigint_cont = 0;
}
#endif

void tv_denoise_core(double *x,double *y,double delta,double eps,double L,double mu,int m,int n,int maxiter,double *kf,double *epsilon_kf){
    
  register double *df,*yk,*wk,*zk,*t,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k, mn = m*n, one = 1;
  double mL = -L;
  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
	register int i1,i2,i3;

  INITBREAK

  df = malloc(m*n*sizeof(double));
  yk = malloc(m*n*sizeof(double));
  wk = malloc(m*n*sizeof(double));
  zk = malloc(m*n*sizeof(double));
  t  = malloc(m*n*sizeof(double));
  uij = malloc(2*sizeof(double));
  
  for (i=0; i<mn; i++){
		wk[i] = 0.0;
		x[i] = y[i];
	}

  for (k=0; k<maxiter; k++) {

    /* step 1 */
    for (i=0; i<mn; i++) df[i] = 0.0;

    pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = x[i1]-x[i3];
				uij[1] = x[i2]-x[i3];
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }

    /*dobj = -s*nrm2(df) + dot(df,y) */
    dobj = ddot_(&mn, df, &one, y, &one) - delta*dnrm2_(&mn, df, &one);

    if (pobj - dobj < eps || STOPMARK) 
      goto cleanup;
    
    /* step 2 */
    /* t = df - L*(x_k-y)
			 y_k = y - min(1/L, s/nrm2(t))*t */

    for (i=0; i<mn; i++) t[i] = df[i];
    daxpy_(&mn, &mL, x, &one, t, &one);
    daxpy_(&mn, &L,  y, &one, t, &one);
    c1 = -1/MAX(L, dnrm2_(&mn, t, &one)/delta);

    for (i=0; i<mn; i++) yk[i] = y[i];
    daxpy_(&mn, &c1, t, &one, yk, &one);

    /* step 3 */
    /* w_k += (k+1)/2.0*df
       z_k = y - min([1/L, s/nrm2(w_k)])*(w_k) */

    c1 = (1.0+k)*(1.0+k)/2.0;
    daxpy_(&mn, &c1, df, &one, wk, &one);
    
    /* for (i=0; i<mn; i++) t[i] = wk[i];
       daxpy_(&mn, &L, y, &one, t, &one); */
    
    for (i=0; i<mn; i++) zk[i] = y[i];
    c1 = -1/MAX(L, dnrm2_(&mn, wk, &one)/delta);
    daxpy_(&mn, &c1, wk, &one, zk, &one);

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2)*(k+2)/2.0;
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;
		
		for (i=0; i<mn; i++) x[i] = t_k*zk[i];
    
		m1t_k = 1-t_k;
    daxpy_(&mn, &m1t_k, yk, &one, x, &one);

  }
  
 cleanup:
  free(df);
  free(yk);
  free(wk);
  free(zk);
  free(t);
  
  kf[0] = (double)(k);
  epsilon_kf[0] = pobj-dobj;

}

void tv_denoise_core_org(double *x,double *y,double delta,double eps,double L,double mu,int m,int n,int maxiter,double *kf,double *epsilon_kf){
//X,B,delta,epsilon,Lmu,mu,size(B,1),size(B,2,),N_maxiter,cur_iter_k,epsilon_k
	
  register double *df,*yk,*wk,*zk,*t,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k, mn = m*n, one = 1;
  double mL = -L;
  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
	register int i1,i2,i3;

  INITBREAK
  
  df = malloc(m*n*sizeof(double));
  yk = malloc(m*n*sizeof(double));
  wk = malloc(m*n*sizeof(double));
  zk = malloc(m*n*sizeof(double));
  t  = malloc(m*n*sizeof(double));
  uij = malloc(2*sizeof(double));
  
  for (i=0; i<mn; i++){
		wk[i] = 0.0;
		x[i] = y[i];//X = B;
	}

  for (k=0; k<maxiter; k++) {

    /* step 1 */
    for (i=0; i<mn; i++) df[i] = 0.0;

    pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = x[i1]-x[i3]; //(Xc')_ij
				uij[1] = x[i2]-x[i3];//(Xr')_ij
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);//TV_mu
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }

    /*dobj = -s*nrm2(df) + dot(df,y) */
    dobj = ddot_(&mn, df, &one, y, &one) - delta*dnrm2_(&mn, df, &one);

    if (pobj - dobj < eps || STOPMARK) 
      goto cleanup;
    
    /* step 2 */
    /* t = df - L*(x_k-y)
			 y_k = y - min(1/L, s/nrm2(t))*t */

    for (i=0; i<mn; i++) t[i] = df[i];
    daxpy_(&mn, &mL, x, &one, t, &one);
    daxpy_(&mn, &L,  y, &one, t, &one);
    c1 = -1/MAX(L, dnrm2_(&mn, t, &one)/delta);

    for (i=0; i<mn; i++) yk[i] = y[i];
    daxpy_(&mn, &c1, t, &one, yk, &one);

    /* step 3 */
    /* w_k += (k+1)/2.0*df
       z_k = y - min([1/L, s/nrm2(w_k)])*(w_k) */

    c1 = (1.0+k)/2.0;
    daxpy_(&mn, &c1, df, &one, wk, &one);
    
    /* for (i=0; i<mn; i++) t[i] = wk[i];
       daxpy_(&mn, &L, y, &one, t, &one); */
    
    for (i=0; i<mn; i++) zk[i] = y[i];
    c1 = -1/MAX(L, dnrm2_(&mn, wk, &one)/delta);
    daxpy_(&mn, &c1, wk, &one, zk, &one);

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2)/2.0;
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;
		
		for (i=0; i<mn; i++) x[i] = t_k*zk[i];
    
		m1t_k = 1-t_k;
    daxpy_(&mn, &m1t_k, yk, &one, x, &one);

  }
  
 cleanup:
  free(df);
  free(yk);
  free(wk);
  free(zk);
  free(t);
  
  kf[0] = (double)(k);
  epsilon_kf[0] = pobj-dobj;

}

void tv_inpaint_core(double *x,double *y,int *Ii,int *Ic,double gamma,double d,double eps,double L,double mu,int m,int n,int sI,int sIc,int maxiter,double *kf, double *epsilon_kf){

  register double *df,*yk,*wk,*zk,*t1,*t2,*yic,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k, mn = m*n, one = 1;
  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
  double Ld = L*d;
	register int i1,i2,i3;

  INITBREAK

  df = malloc(m*n*sizeof(double));
  yk = malloc(sI*sizeof(double));
  wk = malloc(sI*sizeof(double));
  zk = malloc(sI*sizeof(double));
  t1  = malloc(sI*sizeof(double));
  t2  = malloc(sIc*sizeof(double));
  yic = malloc(sIc*sizeof(double));
  uij = malloc(2*sizeof(double));

  for (i=0; i<mn; i++){
		x[i] = y[i];
	}

  for (i=0; i<sIc; i++){
		yic[i] = y[Ic[i]];
	}
	
	for (i=0; i<sI; i++){
		wk[i] = 0.0;
	}

  for (k=0; k<maxiter; k++) {

    /* step 1 */
    for (i=0; i<mn; i++) df[i] = 0.0;
    
    pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = x[i1]-x[i3];
				uij[1] = x[i2]-x[i3];
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }
    
    /* #dobj = -gamma*nrm2(I*D.T*u) + (u.T*D*Ic.T*Ic*b)+ (u.T*D*I.T*d)
			 dobj = -gamma*nrm2(df[I]) + dot(df[Ic],b[Ic]) +dot(df[I],d)*/

		dobj = 0;
		for (i=0; i<sI; i++){			
			t1[i] = df[Ii[i]];
			dobj += t1[i];
		}
		

		dobj=dobj*d;

		for (i=0; i<sIc; i++)
			t2[i] = df[Ic[i]];

    dobj += ddot_(&sIc, t2, &one, yic, &one) - gamma*dnrm2_(&sI, t1, &one);

    if (pobj - dobj < eps || STOPMARK) 
      goto cleanup;
    
    /* step 2 */
		/* t1 = L*x_k[I] - df[I] -L*d
       y_k[I]  = min(1/L, gamma/nrm2(t1))*t1 + d*/

    for (i=0; i<sI; i++) 			
			t1[i] = L*x[Ii[i]] - df[Ii[i]]-Ld;

    c1 = 1/MAX(L, dnrm2_(&sI, t1, &one)/gamma);

		for (i=0; i<sI; i++)
			yk[i] = c1*t1[i]+d;


    /* step 3 */
    /* w_k += (k+1)/2.0*df_I */
    c1 = (1.0+k)/2.0;

		for (i=0; i<sI; i++)
			t1[i] = df[Ii[i]];

    daxpy_(&sI, &c1, t1, &one, wk, &one);
		
		/* t1 = -w_k[I]
       z_k[I]  = min(1/L, gamma/nrm2(t1))*t1+d*/
		
		for (i=0; i<sI; i++) 			
			t1[i] = -wk[i];

    c1 = 1/MAX(L, dnrm2_(&sI, t1, &one)/gamma);

		for (i=0; i<sI; i++) 			
			zk[i] = c1*t1[i]+d;

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2)/2.0;
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;
		
		for (i=0; i<sI; i++) t1[i] = t_k*zk[i];
    
		m1t_k = 1-t_k;
    daxpy_(&sI, &m1t_k, yk, &one, t1, &one);
		
		/*Update the inpainted pixels*/
		for (i=0; i<sI; i++) x[Ii[i]] = t1[i];

  }
  
 cleanup:
  free(df);
  free(yk);
  free(wk);
  free(zk);
  free(t1);
  free(t2);
  free(yic);

  kf[0] = (double)(k);
  epsilon_kf[0] = pobj-dobj;

}

void tv_denoise_inpaint_core(double *x,double *y,int *Ii,int *Ic,double delta,double gamma,double d,double eps,double L,double mu,int m,int n,int sI,int sIc,int maxiter,double *kf, double *epsilon_kf){

  register double *df,*yk,*wk,*zk,*t1,*t2,*yic,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k, mn = m*n, one = 1;
  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
  double Ld = L*d;
	register int i1,i2,i3;

  INITBREAK

  df = malloc(m*n*sizeof(double));
  yk = malloc(m*n*sizeof(double));
  wk = malloc(m*n*sizeof(double));
  zk = malloc(m*n*sizeof(double));
  t1  = malloc(sI*sizeof(double));
  t2  = malloc(sIc*sizeof(double));
  yic = malloc(sIc*sizeof(double));
  uij = malloc(2*sizeof(double));

  for (i=0; i<mn; i++){
		wk[i] = 0.0;
		x[i] = y[i];
	}

  for (i=0; i<sIc; i++)
		yic[i] = y[Ic[i]];

  for (k=0; k<maxiter; k++) {

    /* step 1 */
    for (i=0; i<mn; i++) df[i] = 0.0;
    
		pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = x[i1]-x[i3];
				uij[1] = x[i2]-x[i3];
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }

    
    /* #dobj = -gamma*nrm2(I*D.T*u)-delta*nrm2(Ic*D.T*u) + (u.T*D*Ic.T*Ic*b)+ (u.T*D*I.T*d)
			 dobj = -gamma*nrm2(df[I]) -delta*nrm2(df[Ic]) + dot(df[Ic],b[Ic]) +dot(df[I],d)*/

		dobj = 0;
		for (i=0; i<sI; i++){			
			t1[i] = df[Ii[i]];
			dobj += t1[i];
		}
		
		dobj=dobj*d;

    for (i=0; i<sIc; i++)
			t2[i] = df[Ic[i]];

    dobj += ddot_(&sIc, t2, &one, yic, &one) - gamma*dnrm2_(&sI, t1, &one) - delta*dnrm2_(&sIc, t2, &one);

    if (pobj - dobj < eps || STOPMARK) 
      goto cleanup;
    
    /* step 2 */
		/* t1 = L*x_k[I] - df[I] -L*d
       t2 = L*(x_k[Ic]-b[Ic]) - df[Ic]
       y_k[I]  = min(1/L, gamma/nrm2(t1))*t1 + d
       y_k[Ic] = min(1/L, delta/nrm2(t2))*t2 + b[Ic] */

    for (i=0; i<sI; i++) 			
			t1[i] = L*x[Ii[i]] - df[Ii[i]]-Ld;

    for (i=0; i<sIc; i++)
			t2[i] = L*(x[Ic[i]]-y[Ic[i]]) - df[Ic[i]];

    c1 = 1/MAX(L, dnrm2_(&sI, t1, &one)/gamma);
		c2 = 1/MAX(L, dnrm2_(&sIc, t2, &one)/delta);

		for (i=0; i<sI; i++)
			yk[Ii[i]] = c1*t1[i]+d;

    for (i=0; i<sIc; i++)
			yk[Ic[i]] = c2*t2[i]+y[Ic[i]];

    /* step 3 */
    /* w_k += (k+1)/2.0*df */
    c1 = (1.0+k)/2.0;
    daxpy_(&mn, &c1, df, &one, wk, &one);
		
		/* t1 = -w_k[I]
       t2 = -w_k[Ic] 
       z_k[I]  = min(1/L, gamma/nrm2(t1))*t1+d
       z_k[Ic] = b[Ic] + min(1/L, delta/nrm2(t2))*t2 */
		
		for (i=0; i<sI; i++) 			
			t1[i] = -wk[Ii[i]];

    for (i=0; i<sIc; i++)
			t2[i] = -wk[Ic[i]];

    c1 = 1/MAX(L, dnrm2_(&sI, t1, &one)/gamma);
		c2 = 1/MAX(L, dnrm2_(&sIc, t2, &one)/delta);

		for (i=0; i<sI; i++) 			
			zk[Ii[i]] = c1*t1[i]+d;

    for (i=0; i<sIc; i++)
			zk[Ic[i]] = y[Ic[i]]+c2*t2[i];

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2)/2.0;
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;
		
		for (i=0; i<mn; i++) x[i] = t_k*zk[i];
    
		m1t_k = 1-t_k;
    daxpy_(&mn, &m1t_k, yk, &one, x, &one);

  }
  
 cleanup:
  free(df);
  free(yk);
  free(wk);
  free(zk);
  free(t1);
  free(t2);
  free(yic);

  kf[0] = (double)(k);
  epsilon_kf[0] = pobj-dobj;

}

#ifdef FFTW3

void tv_deblur_core(double *x,double *y,double *s,double delta,double eps, double L, double mu,int m,int n,int maxiter,double *kf, double *epsilon_kf){

  register double *df,*yk,*wk,*zk,*t,*t2,*q,*qs,*dis,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k=0, mn = m*n, one = 1;
  double miL = - 1/L;

  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
  double nuc, *nu1,*nu2;
  double nu1v=0,nu2v=0,mins=1e16,nrm2y;
	register int i1,i2,i3;

  fftw_plan *p1,*p2;
  
  INITBREAK

  df = malloc(m*n*sizeof(double));
  yk = malloc(m*n*sizeof(double));
  wk = malloc(m*n*sizeof(double));
  zk = malloc(m*n*sizeof(double));
  t  = malloc(m*n*sizeof(double));
  t2 = malloc(m*n*sizeof(double));
  q  = malloc(m*n*sizeof(double));
  qs = malloc(m*n*sizeof(double));
  dis = malloc(m*n*sizeof(double));
  uij = malloc(2*sizeof(double));
  
  p1 = malloc(sizeof(fftw_plan));
  p2 = malloc(sizeof(fftw_plan));

  nu1 = &nu1v;
  nu2 = &nu2v;
	
  odct2(y,m,n,p1,0); /* Inplace */

  for (i=0; i<mn; i++){
		wk[i] = 0.0;
		x[i] = y[i];

		/* also set the dis used for trp in step 3 and 4 */
    dis[i] = 1/(s[i]*s[i]);

	}

	nrm2y = dnrm2_(&mn,y,&one);

	mins = fabs(s[0]);
	for (i=0; i<mn; i++){
		mins = MIN(mins,fabs(s[i]));
	}


	nuc = -1*minf(dis,mn); /* max(-x) = -min(x) */


  for (k=0; k<maxiter; k++) {

		
    /* step 1 */
    for (i=0; i<mn; i++){
			df[i] = 0.0;
			t[i] = x[i];
		}

		oidct2(t,m,n,p2,k); /* Inplace oidct2 transform */

		pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = t[i1]-t[i3];
				uij[1] = t[i2]-t[i3];
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }

		odct2(df,m,n,p1,k); /* Inplace */
	
    /* dobj = -delta*nrm2(u'*A*C'*Di) + u.T*A*C'*Di*yt
      			 = -delta*nrm2(Di*df) + yt'*Di*df
						 dobj = -delta*nrm2(div(df,S)) + dot(Y,div(df,S)); */

		for (i=0; i<mn; i++){
			t[i] = df[i]/s[i];
		}

		dobj = ddot_(&mn, t, &one, y, &one) - delta*dnrm2_(&mn, t, &one);

		if (pobj - dobj < eps || STOPMARK)
			goto cleanup;
    
    /* step 2 */
		/* T = X_k - df/L
			y_k nu1 = trp(t, s, y, delta, nu1)
		*/
    for (i=0; i<mn; i++) t[i] = x[i];

    daxpy_(&mn, &miL, df, &one, t, &one);

		trp(t,s,dis,nuc,y,delta,nu1,t2,qs,q,mn,nrm2y,mins,yk); /* solution ends in yk */

    /* step 3 */
    /* w_k += (k+1)/2.0*df */

    c1 = (1.0+k)/2.0;
    daxpy_(&mn, &c1, df, &one, wk, &one);
    
		/*z_k, nu2 = trp(- w_k/L, s, y, delta, nu2)*/
	
    for (i=0; i<mn; i++) t[i] = miL*wk[i];

		trp(t,s,dis,nuc,y,delta,nu2,t2,qs,q,mn,nrm2y,mins,zk);

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2) /2.0; /* alternatively a power 1.19 */
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;
		
		for (i=0; i<mn; i++) x[i] = t_k*zk[i];
    
		m1t_k = 1-t_k;
    daxpy_(&mn, &m1t_k, yk, &one, x, &one);

		/*if(k==2)
			goto cleanup; */
  }
  
 cleanup:
  oidct2(x,m,n,p2,0);
  oidct2(y,m,n,p2,0);

  free(df);
  free(yk);
  free(wk);
  free(zk);
  free(t);
  free(t2);
  free(q);
  free(qs); 
  free(dis);
  fftw_destroy_plan(p1[0]);
  fftw_destroy_plan(p2[0]);
    
  kf[0] = (double)(k);
  epsilon_kf[0] = pobj-dobj;
}



void tv_deblur_rr_core(double *x,double *y,double *s,int *J,int *Jc,double delta,double gamma,double eps, double L, double mu,int m,int n,int sJ,int sJc,int maxiter,double *kf, double *epsilon_kf){

  register double *df,*dfJ,*dfJc,*ykJ,*ykJc,*wk,*zkJ,*zkJc,*t,*tJ,*tJc,*tJ2,*spJ,*yJ,*qJ,*qsJ,*disJ,*uij;

  double pobj = 0, dobj = 0, c1, c2;
  int i, j, k=0, mn = m*n, one = 1;
  double miL = - 1/L,mL=-L,iL=1.0/L;

  double A_kp1=0.5,alpha_kp1,t_k,m1t_k;
  double nuc, *nu1,*nu2;
  double nu1v=0,nu2v=0,minJs=1e16,nrm2Jy=0.0;
	register int i1,i2,i3;

  fftw_plan *p1,*p2;
  
  INITBREAK
  
  df = malloc(m*n*sizeof(double));
  dfJ = malloc(sJ*sizeof(double));
  dfJc = malloc(sJc*sizeof(double));

  ykJ = malloc(sJ*sizeof(double));
  ykJc = malloc(sJc*sizeof(double));

  wk = malloc(m*n*sizeof(double));
  zkJ = malloc(sJ*sizeof(double));
  zkJc = malloc(sJc*sizeof(double));

  t  = malloc(m*n*sizeof(double));
  tJ  = malloc(sJ*sizeof(double));
  tJc = malloc(sJc*sizeof(double));
  tJ2 = malloc(sJ*sizeof(double));

  spJ = malloc(sJ*sizeof(double));
  yJ = malloc(sJ*sizeof(double));

  qJ  = malloc(sJ*sizeof(double));
  qsJ = malloc(sJ*sizeof(double));
  disJ = malloc(sJ*sizeof(double));
  
  uij = malloc(2*sizeof(double));
  
  p1 = malloc(sizeof(fftw_plan));
  p2 = malloc(sizeof(fftw_plan));

  nu1 = &nu1v;
  nu2 = &nu2v;
	 
  odct2(y,m,n,p1,0); /* Inplace */

	for (i=0; i<sJ; i++){
			spJ[i] = s[J[i]];
			yJ[i] = y[J[i]];
			/* also set the dis used for trp in step 3 and 4 */
			disJ[i] = 1/(s[J[i]]*s[J[i]]);
	}
	
  /*calculte some parameters used to determine TOL in trp */

	nrm2Jy = dnrm2_(&sJ,yJ,&one);
	
	minJs = fabs(spJ[0]);
	for (i=1; i<sJ; i++){
		minJs = MIN(minJs,fabs(spJ[i]));
	}

	/* make sure w_k = 0 */
  for (i=0; i<mn; i++){
		wk[i] = 0.0;
		x[i] = y[i];
	}

	nuc = -1.0*minf(disJ,sJ); /* max(-x) = -min(x) */

  for (k=0; k<maxiter; k++) {

    /* step 1 */
    for (i=0; i<mn; i++){
			df[i] = 0.0;
			t[i] = x[i];
		}

		oidct2(t,m,n,p2,k); /* Inplace oidct2 transform */

    pobj = 0.0;
    for (j=0; j<n-1; j++) {
      for (i=0; i<m-1; i++) {
				i1 = (i+1) + j*m;
				i2 = i + (j+1)*m;
				i3= i+j*m;

				uij[0] = t[i1]-t[i3];
				uij[1] = t[i2]-t[i3];
			
				c1 = sqrt(uij[0]*uij[0] + uij[1]*uij[1]);
				pobj += c1; 

				c2 = MAX(mu, c1);
				uij[0] = uij[0]/c2;
				uij[1] = uij[1]/c2;

				df[i1] += uij[0];
				df[i3] -= uij[0];
				df[i2] += uij[1];
				df[i3] -= uij[1];
      }
    }

		odct2(df,m,n,p1,k); /* Inplace */

    /*  dobj = -deltad*nrm2(div(df[I],S[I]))-gamma*nrm2(df[Ic]) + dot(B[I],div(df[I],S[I])) */

		for (i=0; i<sJ; i++){
			tJ[i] = df[J[i]]/spJ[i];
		}

		for (i=0; i<sJc; i++){
			tJc[i] = df[Jc[i]];
		}

		dobj = ddot_(&sJ, tJ, &one, yJ, &one) - delta*dnrm2_(&sJ, tJ, &one)- gamma*dnrm2_(&sJc,tJc,&one);

		if (pobj - dobj < eps || STOPMARK)			
			goto cleanup;
        
    /* step 2 
		 T = X_k_J - df_J/L
			y_k nu1 = trp(t, s, y, delta, nu1)
			#For the complementary part
			t = L*X_k[Jc]-df[Jc]
			Y_k[Jc] = min(1/L, gamma/nrm2(t))*t
		*/

    for (i=0; i<sJ; i++){
			tJ[i] = x[J[i]];
			dfJ[i] = df[J[i]];
		}

    daxpy_(&sJ, &miL, dfJ, &one, tJ, &one);

		trp(tJ,spJ,disJ,nuc,yJ,delta,nu1,tJ2,qsJ,qJ,sJ,nrm2Jy,minJs,ykJ); /* solution ends in yk */

		for (i=0; i<sJc; i++){
			tJc[i] = x[Jc[i]];
			ykJc[i] = -df[Jc[i]];
		}

		daxpy_(&sJc,&L,tJc,&one,ykJc,&one);

		c1 = 1/MAX(L, dnrm2_(&sJc, ykJc, &one)/gamma);
    for (i=0; i<sJc; i++) ykJc[i] = c1*ykJc[i];

    /* step 3 */
    /* w_k += (k+1)/2.0*df */

    c1 = (1.0+k)/2.0;

    daxpy_(&mn, &c1, df, &one, wk, &one);
    
		/*z_k, nu2 = trp(- w_k_J/L, s_J, y_J, delta, nu2)*/

		/*  For the complementay part
				t = - W_k[Ic]
        Z_k[Ic] = min(1/L, gamma/nrm2(t))*t */

    for (i=0; i<sJ; i++) tJ[i] = miL*wk[J[i]];

		trp(tJ,spJ,disJ,nuc,yJ,delta,nu2,tJ2,qsJ,qJ,sJ,nrm2Jy,minJs,zkJ);

		for (i=0; i<sJc; i++) tJc[i] = -wk[Jc[i]];
		c1 = 1/MAX(L, dnrm2_(&sJc, tJc, &one)/gamma);
		for (i=0; i<sJc; i++) zkJc[i] = c1*tJc[i];

    /* step 4 */
    /* x_k = z_k*2/(k+3) + y_k*(k+1)/(k+3) */
		alpha_kp1 = (k+2) /2.0; 
		A_kp1 += alpha_kp1;

		t_k = alpha_kp1/A_kp1;

		m1t_k = 1-t_k;

		for (i=0; i<sJ; i++) x[J[i]] = t_k*zkJ[i] + m1t_k*ykJ[i];
 		for (i=0; i<sJc; i++) x[Jc[i]] = t_k*zkJc[i] + m1t_k*ykJc[i];

		for (i=0; i<sJc; i++)
			tJc[i] =  x[Jc[i]];


		/* Check if gamma bound is reached, then quit */

			 if(gamma < 1.001*dnrm2_(&sJc,tJc,&one))
				 goto cleanup;
  }

 cleanup:
	
  oidct2(x,m,n,p2,0);
  oidct2(y,m,n,p2,0);

  free(df);free(dfJ);free(dfJc);
  free(qJ);	free(qsJ);	free(disJ);
  free(uij);

  free(ykJ); free(ykJc);
  free(wk);
  free(zkJ); free(zkJc);

  free(t);free(tJ);  free(tJc);  free(tJ2);
  free(spJ);
  free(yJ);
	
  fftw_destroy_plan(p1[0]);fftw_destroy_plan(p2[0]);
	
  kf[0] = (double) k;
  epsilon_kf[0] = pobj-dobj;
}

#endif

