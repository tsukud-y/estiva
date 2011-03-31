#include <stdio.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"

static void cpy(double *src, double *dst)
{
  long i;
  forall (0, i, dim1(dst) ) dst[i] = src[i];
}

static int addvec(double da, double *dx,  double *dy)
{
  long i;
  forall (0, i, dim1(dy) ) dy[i] += da * dx[i];
  return 0;
}

static double dotvec(double *dx, double *dy)
{
  double tmp = 0.0;
  long i;
  forall (0, i, dim1(dy) ) tmp += dy[i] * dx[i];  
  return tmp;
}

static double L2(double *dx)
{
  double sum = 0.0;
  long i;
  forall (0, i, dim1(dx) ) sum += dx[i]*dx[i];
  return sqrt(sum);
}

static void matvec(MX *A, double alpha, double *p, double beta, double *q)
{
  long i, j, m, n, J;
  mx(A,1,1) = mx(A,1,1);

  m = A->m;
  n = A->n;
  forall (0, i, m) q[i] = 0.0;

  forall (0, i, m-1) forall(0, j, n-1) {
      J = A->IA[i][j];
      if (J != 0) q[i+1] += A->A[i][j]*p[J];
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

static int success(long k)
{
  if ( defop("-v") ) printf("itr = %ld\n",k);
  return 0;
}

static void phase(int i)
{
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

static void ILU(MX *A, double *d, double *dd, double *b)
{
  double ss;
  long   i, k, n, K, An;

  mx(A,1,1) = mx(A,1,1);

  An = A->n;
  n  = dim1(b);

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

int estiva_pcgssolver(void *A, double *xk, double *b)
{
  static double *d, *dd, *ek, *ek1, *hk1, *pk, *pk1, *q, *r0, *rk, *rk1, 
    *w, *xk1, *xk1_xk;
  double c1, c2, c3, alphak, betak, epsilon = 1.0e-7;
  long   i, k, n;

  n = dim1(b);

  ary1(d  ,n+1);
  ary1(dd ,n+1);
  ary1(ek ,n+1);
  ary1(ek1,n+1);
  ary1(hk1,n+1);
  ary1(pk ,n+1);
  ary1(pk1,n+1);
  ary1(q  ,n+1);
  ary1(r0 ,n+1);
  ary1(rk ,n+1);
  ary1(rk1,n+1);
  ary1(w  ,n+1);
  ary1(xk1,n+1);
  ary1(xk1_xk,n+1);

  forall (1, i, n) d[i] = mx(A,i,i);   
  ILU(A,d,dd,b);
  matvec(A,1.0,xk,0.0,q);
  formula( rk, '=', b, '-', 1.0,q);
  presolve(A,dd,rk);
  cpy(rk,r0);
  cpy(rk,ek);
  cpy(rk,pk);
  c1 = dotvec(r0,r0);

  forall (1, k, n) {
    phase(1); {
      matvec(A, 1.0, pk, 0.0, q);
      presolve(A,dd,q);
      c2 = dotvec(q,r0);
      alphak = c1 / c2;
    }
    phase(2); {
      formula( hk1, '=', ek,'-',alphak,q);
      formula( w,   '=', ek,'+',1.0,hk1); 
    }
    phase(3); {
      matvec(A, 1.0, w, 0.0, q);
      presolve(A,dd,q);
    }
    phase(4); {
      formula( rk1, '=', rk,'-',alphak,q);
      formula( xk1, '=', xk,'+',alphak,w);
      c3 = dotvec(rk1,r0);
    }
    phase(5); {
      formula( xk1_xk, '=', xk1,'-',1.0,xk);
      if ( L2(xk1_xk)/L2(xk) < epsilon ) return success(k);
    }
    phase(6); {
      betak = c3 / c1;
      c1 = c3;
      formula( ek1, '=', rk1,'+',betak,hk1);
      formula( pk1, '=', ek1,'+',betak, formula( w, '=', hk1,'+',betak,pk));
    }
    cpy(ek1,ek);
    cpy(pk1,pk);
    cpy(rk1,rk);
    cpy(xk1,xk);
  }
  return 1;
}
