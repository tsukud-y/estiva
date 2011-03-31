#include <math.h>
#include <stdio.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"

static void ILUdecomposition(MX *A, double *dd, double *b)
{
  static double *d;
  double ss;
  long   i, k, n, K, An;
  An = A->n;
  n  = dim1(b);
  ary1(d,n+1);

  forall (1, i, n) d[i] = mx(A,i,i);

  mx(A,1,1) = mx(A,1,1);

  dd[1] = 1.0 / d[1];
  forall (2, i, n) {
    ss = d[i];
    forall ( 0, k, An-1) {
      K = A->IA[i-1][k];
      if ( K != 0 && i > K ) {
        ss -= A->A[i-1][k] * mx(A,K,i) * dd[K];
      }
    }
    dd[i] = 1.0 / ss;
  }
}

static void presolve(MX *A, double *dd, double *q)
{
  double sw;
  long   i, j, n, J, An;

  mx(A,1,1) = mx(A,1,1);
  n  = dim1(q);
  An = A->n;

  forall (1, i, n) {
    forall (0, j, An-1) {
      J = A->IA[i-1][j];
      if ( J != 0 && i > J ) 
        q[i] -= A->A[i-1][j]*q[J];
    }
    q[i] *= dd[i];
  }
  for (i=n; i>0; i--) {
    sw = 0.0;
    forall (0, j, An-1) {
      J = A->IA[i-1][j];
      if ( J != 0 && i < J )
        sw += A->A[i-1][j]*q[J];
    }
    q[i] -= dd[i] * sw;
  }
}

static int matvec(MX *A, double alpha, double *x, double beta, double *y)
{ 
  matvecmx(A, &alpha, &x[1], &beta, &y[1]); 
  return 0;
}

static int addvec(double da, double *dx,  double *dy)
{
  long i ;
  forall (0, i, dim1(dy) ) dy[i] += da * dx[i];
  return 0;
}

static double dotvec(double *dx, double *dy)
{
  double sum = 0.0;
  long i;
  forall (0, i, dim1(dy) ) sum += dy[i] * dx[i];  
  return sum;
}

static double L2(double *dx)
{
  double sum = 0.0;
  long i;
  forall (0, i, dim1(dx) ) sum += dx[i]*dx[i];
  return sqrt(sum);
}

static void cpy(double *src, double *dst)
{
  long i;
  forall (0, i, dim1(dst) ) dst[i] = src[i];
}

static double *formula(double *z, char eq, 
                       double *x, char plus, double a, double *y)
{
  cpy(x,z);
  if ( plus == '+' )
    addvec(a,y,z);
  else
    addvec(-a,y,z);
  return z;
}

int estiva_cgssolver(void *A, double *x, double *b)
{
  static double *dd, *p, *phat, *q, *qhat, *r, *rtld, *tmp, *u, *uhat, *vhat;
  double alpha, beta, bnrm2, rho, rho1=1.0;
  long   itr, n;

  n = dim1(b);

  ary1(r    ,n+1);
  ary1(rtld ,n+1);
  ary1(p    ,n+1);
  ary1(phat ,n+1);
  ary1(q    ,n+1);
  ary1(qhat ,n+1);
  ary1(u    ,n+1);
  ary1(uhat ,n+1);
  ary1(vhat ,n+1);
  ary1(tmp  ,n+1);
  ary1(dd,   n+1);

  cpy(b,x);
  ILUdecomposition(A,dd,b);
  cpy(b, r);
  if ( L2(x) != 0.0 ) {
    matvec(A,-1.0, x, 1.0, r);
    if (L2(r) <= 1.0e-7) return 0;
  }
  bnrm2 = L2(b);
  if (bnrm2 == 0.0) bnrm2 = 1.0;
  cpy(r, rtld);

  for (itr = 1; itr < n;  itr++) {
    rho = dotvec(rtld,r);
    if (fabsl(rho) < 1.2e-31) break;
    if ( itr == 1 ) {
      cpy(r,u);
      cpy(u,p);
    } else {
      beta = rho / rho1;
      formula( u, '=', r, '+', beta, q);
      formula( p, '=', u, '+', beta, formula( tmp, '=',q,'+',beta,p));
    }
    cpy(p,phat);
    presolve(A,dd,phat);
    matvec(A,1.0, phat, 0.0, vhat );

    alpha = rho / dotvec(rtld, vhat);
    formula( q, '=', u, '-', alpha, vhat);
    cpy(formula(phat, '=', u, '+', 1.0, q),uhat);
    presolve(A,dd,uhat);

    formula(x, '=', x, '+', alpha, uhat);

    matvec(A,1.0, uhat, 0.0, qhat );
    formula(r, '=', r, '-',alpha,qhat);

    if ( L2(r) / bnrm2 <= 1.0e-7) break;
    rho1 = rho;
  }
  if ( defop("-v") ) printf("itr = %ld\n",itr);
  return 0;
}
