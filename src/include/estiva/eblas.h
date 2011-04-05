#ifndef _ESTIVA_EBLAS_H_
#define _ESTIVA_EBLAS_H_

int    estiva_daxpy_(int *n, double *alpha, double *x, int *c1, double *y, int *c2);
int    estiva_dcopy_(int *n,double  *b, int *c__1, double *r, int *c__2);
double estiva_ddot_(int *n, double *z, int *c__1, double *y, int *c__2);
double estiva_dnrm2_(int *n, double *x, int *c__1);
int    estiva_dscal_(int *n, double *d__1, double *v, int *c__1);
double estiva_getbreak_(void);
void  *estiva_setAmx(void *Apointer);
void  *estiva_setATmx(void *ATpointer);
int    estiva_matvec(double *alpha, double *x, double *beta, double *y);
int    estiva_matvectrans(double *alpha, double *x, double *beta, double *y);
int    estiva_psolveq(double *x, double *b, char *str, int L);
int    estiva_psolvetransq(double *x, double *b, char *str, int L);

#define daxpy_(n,alpha,x,c1,y,c2) estiva_daxpy_(n,alpha,x,c1,y,c2)
#define dcopy_(n,b,c__1,r,c__2)   estiva_dcopy_(n,b,c__1,r,c__2)
#define ddot_(n,z,c__1,y,c__2)    estiva_ddot_(n,z,c__1,y,c__2)
#define dnrm2_(n,x,c__1)          estiva_dnrm2_(n,x,c__1)
#define dscal_(n,d__1,v,c__1)     estiva_dscal_(n,d__1,v,c__1)
#define getbreak_()               estiva_getbreak_()
#define setAmx(A)                 estiva_setAmx(A)
#define setATmx(AT)               estiva_setATmx(AT)

#endif
