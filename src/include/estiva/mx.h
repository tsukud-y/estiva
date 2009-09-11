#ifndef _ESTIVA_MX_H_
#define _ESTIVA_MX_H_

typedef struct{
  double **A;
  long   **IA;
  long   m, n, I, J;
  double a;
} MX;

#define initmx(A,m,n)         estiva_initmx(&A,m,n)
#define mx(A,i,j)           (*estiva_mx(A,i,j))
#define fmx(A)                estiva_fmx(A)
#define rmx(A,i,j)            estiva_rmx(A,i,j)
#define wmx(A,i,j,a)          estiva_write_matrix(A,i,j,a)


extern void    estiva_initmx(MX **Ap, long i, long j);
extern double *estiva_mx(MX *T, long i, long j);
extern void    estiva_fmx(MX *T);
extern double  estiva_rmx(MX *T, long i, long j);
extern void    estiva_write_matrix(MX *T, long i, long j, double a);


#endif
