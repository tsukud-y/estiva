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

void    estiva_initmx(MX **Ap, long i, long j);
void    estiva_transmx(MX **ATp, MX *A);
void    estiva_clonemx(MX **ATp, MX *A);
void    estiva_mulmx(double **tp, MX *A, double *x);
double *estiva_mx(MX *T, long i, long j);
void    estiva_matvecmx(MX *A, double *alpha, double *x, double *beta, double *y);
void    estiva_psolvemx(MX *A, CRS *pivot, MX *LU, double *D, double *x, double *b);
void    estiva_pltmx(MX *A);
void    estiva_clearmx(MX *A);
void    estiva_fornonzeromx_init(MX *A);
int     estiva_fornonzeromx_loop(MX *A, long *Ip, long *Jp);
void    estiva_zerofillrow(MX *A, long i);

#define initmx(A,m,n)              estiva_initmx(&A,m,n)
#define transmx(AT,A)              estiva_transmx(&AT,A)
#define clonemx(AT,A)              estiva_clonemx(&AT,A)
#define mulmx(t,A,x)               estiva_mulmx(&t,A,x)
#define mx(A,i,j)                (*estiva_mx(A,i,j))
#define matvecmx(A,alpha,x,beta,y) estiva_matvecmx(A,alpha,x,beta,y)
#define psolvemx(A,pivot,LU,D,x,b) estiva_psolvemx(A,pivot,LU,D,x,b)
#define pltmx(A)                   estiva_pltmx(A)
#define clearmx(A)                 estiva_clearmx(A)
#define fornonzeromx(A,i,j) for(estiva_fornonzeromx_init(A);estiva_fornonzeromx_loop(A,&(i),&(j));)
#define zerofillrow(A,i)           estiva_zerofillrow(A,i)

#endif
