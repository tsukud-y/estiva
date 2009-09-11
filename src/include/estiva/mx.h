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
#define write_matrix(A,i,j,a) estiva_write_matrix(A,i,j,a)
#define read_matrix(A,i,j)    estiva_read_matrix(A,i,j)
#define flush_matrix(A)       estiva_flush_matrix(A)

extern void    estiva_initmx(MX **Ap, long i, long j);
extern double *estiva_mx(MX *T, long i, long j);
extern double  estiva_read_matrix(MX *T, long i, long j);
extern void    estiva_write_matrix(MX *T, long i, long j, double a);
extern void    estiva_flush_matrix(MX *T);

#endif
