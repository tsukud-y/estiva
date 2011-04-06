#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"

static long dim1vec;
static double *DD;


void estiva_setveclength(long n)
{
  dim1vec = n-1;
}

double *estiva_cpvec(double *src, double *dst)
{
  long i;
  forall (0, i, dim1vec ) dst[i] = src[i];
  return dst;
}

double *estiva_addvec(double da, double *dx,  double *dy)
{
  long i;
  forall (0, i, dim1vec ) dy[i] += da * dx[i];
  return dy;
}

double *estiva_addformula(double *z, char eq, 
			  double *x, char plus, double a, double *y)
{
  cpvec(x,z);
  if ( plus == '+' )
    ;
  else
    a = -a;
  addvec(a,y,z);
  return z;
}

double estiva_dotvec(double *dx, double *dy)
{
  double sum = 0.0;
  long i;
  forall (0, i, dim1vec ) sum += dy[i] * dx[i];  
  return sum;
}

double estiva_L2(double *dx)
{
  return sqrt(dotvec(dx,dx));
}


double *estiva_matvecvec(void *Apointer, double alpha, double *p, 
			 double beta, double *q)
{
  MX *A;
  A = Apointer;
  mx(A,1,1) = mx(A,1,1);

  if ( A->n == dim1vec ) {
    matvecmx(A,&alpha,p+1,&beta,q+1);
  } else {
    matvecmx(A,&alpha,p,&beta,q);
  }
  return q;
}

double *estiva_psolvevec(void *Apointer, double *q)
{
  MX *A;
  double sw;
  long   i, j, n, J, An;

  A = Apointer;
  n  = A->n;
  An = A->w;

  mx(A,1,1) = mx(A,1,1);

  if ( n == dim1vec ) ;
  else q--;
  
  forall (1, i, n) {
    forall (0, j, An-1) {
      J = A->IA[i-1][j];
      if ( i > J && J > 0) q[i] -= A->A[i-1][j]*q[J];
    }
    q[i] *= DD[i];
  }
  for (i=n; i>0; i--) {
    sw = 0.0;
    forall (0, j, An-1) {
      J = A->IA[i-1][j];
      if ( i < J && 0 < J ) sw += A->A[i-1][j]*q[J];
    }
    q[i] -= DD[i] * sw;
  }
  return q+1;
}

double estiva_Linf(double *x)
{
  double val = 0.0;
  long i;
  
  forall (0, i, dim1vec) {
    val = max(val,fabs(x[i]));
  }
  return val;
}


double *estiva_ILUdecomp(void *Apointer)
{
  static double *d;
  MX *A;
  double ss;
  long   i, k, n, K, An;

  A = Apointer;

  An = A->w;
  n  = A->n;
  
  ary1(d,n+1);
  ary1(DD,n+1);
  
  forall(1,i,n) d[i] = mx(A,i,i);

  mx(A,1,1) = mx(A,1,1);

  DD[1] = 1.0 / d[1];
  forall (2, i, n) {
    ss = d[i];
    forall ( 0, k, An-1) {
      K = A->IA[i-1][k];
      if ( i > K && 0 < K) {
	ss -= A->A[i-1][k] * mx(A,K,i) * DD[K];
      }
    }
    DD[i] = 1.0 / ss;
  }
  return DD;
}

int estiva_success(long k)
{
  if ( defop("-v") ) {
    if (defop("-epsilon")) printf("epsilon = %e\n",epsilon());
    printf("itr = %ld\n",k);
  }
  return 0;
}

double estiva_epsilon(void)
{
  double eps = 1.0e-7;

  if ( defop("-epsilon") ) {
    eps= atof(getop("-epsilon"));
  }
  return eps;
}

int estiva_phase(int i)
{
  return 0;
}

double *estiva_scalvec(double da, double *dx)
{
  long i__2, n;
  static long i, m, mp1;
  n = dim1vec;
  --dx;
  m = n % 5;
  if (m == 0) {
    goto L40;
  }
  i__2 = m;
  for (i = 1; i <= i__2; ++i) {
    dx[i] = da * dx[i];
  }
  if (n < 5) {
    return 0;
  }
 L40:
  mp1 = m + 1;
  i__2 = n;
  for (i = mp1; i <= i__2; i += 5) {
    dx[i] = da * dx[i];
    dx[i + 1] = da * dx[i + 1];
    dx[i + 2] = da * dx[i + 2];
    dx[i + 3] = da * dx[i + 3];
    dx[i + 4] = da * dx[i + 4];
  }
  return dx+1;
}
