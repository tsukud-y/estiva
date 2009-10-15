#ifndef _ESTIVA_MX_H_
#define _ESTIVA_MX_H_

typedef struct{
  double **A;
  long   **IA;
  long   m, n, I, J;
  double a;
} MX;

#define initmx(A,m,n)         estiva_initmx(&A,m,n)
#define transmx(AT,A)         estiva_transmx(&AT,A)
#define mulmx(t,A,x)          estiva_mulmx(&t,A,x)
#define mx(A,i,j)           (*estiva_mx(A,i,j))

extern void    estiva_initmx(MX **Ap, long i, long j);
extern void    estiva_transmx(MX **ATp, MX *A);
extern void    estiva_mulmx(double **tp, MX *A, double *x);
extern double *estiva_mx(MX *T, long i, long j);


#endif
