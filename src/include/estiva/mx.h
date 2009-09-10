#ifndef _ESTIVA_MX_H_
#define _ESTIVA_MX_H_

typedef struct{
  double **A;
  long   **IA;
  long   m, n, I, J;
  double a;
} MATRIX;

typedef MATRIX MX;

#define init_matrix(A,m,n)    estiva_init_matrix(&A,m,n)
#define initmx(A,m,n)         estiva_init_matrix(&A,m,n)
#define write_matrix(A,i,j,a) estiva_write_matrix(A,i,j,a)
#define read_matrix(A,i,j)    estiva_read_matrix(A,i,j)
#define matrix(A,i,j)       (*estiva_matrix(A,i,j))
#define mx(A,i,j)           (*estiva_matrix(A,i,j))
#define flush_matrix(A)       estiva_flush_matrix(A)

extern void    estiva_init_matrix(MATRIX **Ap, long i, long j);
extern double  estiva_read_matrix(MATRIX *T, long i, long j);
extern void    estiva_write_matrix(MATRIX *T, long i, long j, double a);
extern double *estiva_matrix(MATRIX *T, long i, long j);
extern void    estiva_flush_matrix(MATRIX *T);

#endif
