#ifndef _ESTIVA_MX_H_
#define _ESTIVA_MX_H_

typedef struct {
  long   m, n, *col_ind, *row_ptr, *diag_ptr;
  double *val, *pivots;
} CRS;

typedef struct{
  double **A;
  long   **IA;
  long   m, n, I, J;
  double a;
} MX;

void    estiva_clearmx(MX *A);
void    estiva_clonemx(MX **ATp, MX *A);
void    estiva_initmx(MX **Ap, long i, long j);
void    estiva_matvecmx(MX *A,double *alpha,double *x,double *beta,double *y);
void    estiva_mulmx(double **tp, MX *A, double *x);
double *estiva_mx(MX *T, long i, long j);
void    estiva_pltmx(MX *A);
void    estiva_psolvemx(MX *A,CRS *pivot,MX *LU,double *D,double *x,double *b);
void    estiva_transmx(MX **ATp, MX *A);
void    estiva_zerofillrow(MX *A, long i);
void    estiva_slimupmx(MX **Ap, MX *A);
void    estiva_fornonzeromx(MX *A);
int     estiva_fornonzeromx_loop(MX *A, long *Ip, long *Jp);
int     estiva_fprintmx(void *A, char *name);

#define clearmx(A)                   estiva_clearmx(A)
#define clonemx(AT,A)                estiva_clonemx(&AT,A)
#define initmx(A,m,n)                estiva_initmx(&A,m,n)
#define matvecmx(A,alpha,x,beta,y)   estiva_matvecmx(A,alpha,x,beta,y)
#define mulmx(t,A,x)                 estiva_mulmx(&t,A,x)
#define mx(A,i,j)                  (*estiva_mx(A,i,j))
#define pltmx(A)                     estiva_pltmx(A)
#define psolvemx(A,pivot,LU,D,x,b)   estiva_psolvemx(A,pivot,LU,D,x,b)
#define transmx(AT,A)                estiva_transmx(&AT,A)
#define zerofillrow(A,i)             estiva_zerofillrow(A,i)
#define slimupmx(As,A)               estiva_slimupmx(&As,A)
#define fornonzeromx(A,i,j)						\
  for ( estiva_fornonzeromx(A);      estiva_fornonzeromx_loop(A,&(i),&(j));)
#define fprintmx(A)  estiva_fprintmx(A,#A)

#endif
