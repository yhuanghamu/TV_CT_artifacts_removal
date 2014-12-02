#ifndef __TV_CORE__
#define __TV_CORE__

void tv_denoise_core(double *x,double *y,double delta,double eps,double L,double mu,int m,int n,int maxiter,double *kf,double *epsilon_kf);

void tv_denoise_core_org(double *x,double *y,double delta,double eps,double L,double mu,int m,int n,int maxiter,double *kf, double *epsilon_kf);

void tv_denoise_inpaint_core(double *x,double *y,int *J,int *Jc,double delta,double gamma,double d,double eps,double L,double mu,int m,int n,int sJ,int sJc, int maxiter,double *kf, double *epsilon_kf);

void tv_inpaint_core(double *x,double *y,int *J,int *Jc,double gamma,double d,double eps,double L,double mu,int m,int n,int sJ,int sJc, int maxiter,double *kf, double *epsilon_kf);

void tv_deblur_core(double *x,double *y,double *s,double delta,double eps, double L, double mu,int m,int n,int maxiter,double *kf, double *epsilon_kf);

void tv_deblur_rr_core(double *x,double *y,double *s,int *J,int *Jc,double delta,double gamma,double eps, double L, double mu,int m,int n,int sJ,int sJc,int maxiter,double *kf, double *epsilon_kf);


#endif
