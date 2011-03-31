#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

  if (A->m != dim1(b)) { b++; printf("Warrning\n"); abort(); }
  mx(A,1,1) = mx(A,1,1);
  n  = A->m;
  An = A->n;

  ary1(d,n+1);

  forall (1, i, n) d[i] = mx(A,i,i);

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

  if (A->m != dim1(q)) { q++; printf("Warrning\n"); abort(); }
  mx(A,1,1) = mx(A,1,1);
  n  = A->m; 
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

static int matvec(MX *A,double alpha, double *x, double beta, double *y)
{ 
  matvecmx(A, &alpha, &x[1], &beta, &y[1]); 
  return 0;
}

static int addvec(double da, double *dx,  double *dy)
{
  long i;
  forall(0,i,dim1(dy)) dy[i] += da * dx[i];
  return 0;
}

static double dotvec(double *dx, double *dy)
{
  double sum = 0.0;
  long i;
  forall(0,i,dim1(dy)) sum += dy[i] * dx[i];  
  return sum;
}

static double L2(double *dx)
{
  return sqrt(dotvec(dx,dx));
}

static void cpy(double *src, double *dst)
{
  long i;
  forall(0,i,dim1(dst)) dst[i] = src[i];
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


int estiva_bicgsolver(void *Ap, double *x, double *b)
{
  static MX *A, *AT;
  static double *p, *ptld, *q, *qtld, *r, *rtld, *z, *ztld, *dd;
  double alpha, beta, bnrm2, rho, rho1=1.0;
  long   n, itr;
  n = dim1(b);

  ary1(r   ,n+1);
  ary1(rtld,n+1);
  ary1(z   ,n+1);
  ary1(ztld,n+1);
  ary1(p   ,n+1);
  ary1(ptld,n+1);
  ary1(q   ,n+1);    
  ary1(qtld,n+1);
  ary1(dd,  n+1);

  slimupmx(A,Ap);
  transmx(AT,A);
  cpy(b,x);
  
  ILUdecomposition(A,dd,b);
  
  cpy(b,r);
  if (L2(x) != 0.0) {
    matvec(A,-1.0, x, 1.0, r);
    if (L2(r) <= 1.0e-7) {
      return 0;
    }
  }
  cpy(r,rtld);
  bnrm2 = L2(b);
  if (bnrm2 == 0.0) bnrm2 = 1.0;

  forall (1, itr, n) {
    cpy(r,z);
    presolve(A,dd,z);
    cpy(rtld,ztld);
    presolve(AT,dd,ztld);
    rho = dotvec(z, rtld);
    if ( fabsl(rho) < 1.2e-31 ) break;
    if (itr == 1) {
      cpy(z,p);
      cpy(ztld,ptld);
    } else {
      beta = rho / rho1;
      formula(z, '=', z, '+', beta,p );
      cpy(z,p);
      formula(ztld, '=', ztld, '+', beta,ptld );
      cpy(ztld,ptld);
    }
    matvec(A, 1.0, p, 0.0, q);
    matvec(AT,1.0, ptld, 0.0, qtld);
    alpha = rho / dotvec(ptld,q);
    formula(    x, '=', x,    '+', alpha,p    );
    formula(    r, '=', r,    '-', alpha,q    );
    formula( rtld, '=', rtld, '-', alpha,qtld );

    if (L2(r)/bnrm2 <= 1.0e-7) break;

    rho1 = rho;
  }
  if ( defop("-v") ) printf("itr = %ld\n",itr);
  return 0;
}
