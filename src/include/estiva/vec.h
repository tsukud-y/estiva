#ifndef _ESTIVA_VEC_H_
#define _ESTIVA_VEC_H_

double *estiva_addformula(double *z, char eq, double *x, char plus, double a,
			  double *y);
double *estiva_addvec(double da, double *dx, double *dy);
double *estiva_cpvec(double *src, double *dst);
double  estiva_dotvec(double *dx, double *dy);
double  estiva_epsilon(void);
double *estiva_ILUdecomp(void *Apointer);
double  estiva_L2(double *dx);
double  estiva_Linf(double *x);
double *estiva_matvecvec(void *A, double al, double *p,	double be, double *q);
int     estiva_phase(int i);
double *estiva_psolvevec(void *Apointer, double *q);
double *estiva_scalvec(double da,double *dx);
void    estiva_setveclength(long n);
int     estiva_stopcondition(void *A, double *x, double *b);
int     estiva_success(long k);
int     estiva_fprintvec(double *b, char *name);
long    estiva_getveclength(void);

#define addformula(z,eq,x,plus,a,y) estiva_addformula(z,eq,x,plus,a,y)
#define addvec(da,dx,dy)            estiva_addvec(da,dx,dy)
#define cpvec(src,dst)              estiva_cpvec(src,dst)
#define dotvec(dx,dy)               estiva_dotvec(dx,dy)
#define epsilon()                   estiva_epsilon()
#define ILUdecomp(A)                estiva_ILUdecomp(A)
#define L2(dx)                      estiva_L2(dx)
#define Linf(x)                     estiva_Linf(x)
#define matvecvec(A,alpha,p,beta,q) estiva_matvecvec(A,alpha,p,beta,q)
#define phase(i)                    estiva_phase(i)
#define psolvevec(A,q)              estiva_psolvevec(A,q)
#define scalvec(da,dx)              estiva_scalvec(da,dx)
#define setveclength(n)             estiva_setveclength(n)
#define stopcondition(A,x,b)        estiva_stopcondition(A,x,b)
#define success(k)                  estiva_success(k)
#define fprintvec(b)                estiva_fprintvec(b,#b)
#define getveclength()              estiva_getveclength()

#endif
