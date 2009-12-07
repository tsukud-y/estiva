#ifndef _ESTIVA_MX_H_
#define _ESTIVA_MX_H_

typedef struct{
  double **A;
  long   **IA;
  long   m, n, I, J;
  double a;
} MX;

typedef struct {
  long   m, n, *col_ind, *row_ptr, *diag_ptr;
  double *val, *pivots;
} CRS;


#define initmx(A,m,n)         estiva_initmx(&A,m,n)
#define transmx(AT,A)         estiva_transmx(&AT,A)
#define clonemx(AT,A)         estiva_clonemx(&AT,A)
#define mulmx(t,A,x)          estiva_mulmx(&t,A,x)
#define mx(A,i,j)           (*estiva_mx(A,i,j))
#define matvecmx(A,alpha,x,beta,y) estiva_matvecmx(A,alpha,x,beta,y)
#define psolvemx(A,pivot,LU,D,x,b) estiva_psolvemx(A,pivot,LU,D,x,b)

extern void    estiva_initmx(MX **Ap, long i, long j);
extern void    estiva_transmx(MX **ATp, MX *A);
extern void    estiva_clonemx(MX **ATp, MX *A);
extern void    estiva_mulmx(double **tp, MX *A, double *x);
extern double *estiva_mx(MX *T, long i, long j);
extern void    estiva_matvecmx(MX *A, double *alpha, double *x, double *beta, double *y);
extern void    estiva_psolvemx(MX *A, CRS *pivot, MX *LU, double *D, double *x, double *b);

#endif
