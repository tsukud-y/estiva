#include <stdio.h>
#include <string.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/precond.h"
#include "estiva/solver.h"
#include "estiva/std.h"


static MX *A, *AT, *LU, *LUT;
static double *D;
static CRS pivot, pivotT;
static int bicg_();

static int matvec(double alpha, double *x, double beta, double *y)
{ 
  matvecmx(A, &alpha, x, &beta, y); 
  return 0;
}

static int matvectrans(double alpha, double *x, double beta, double *y)
{ 
  matvecmx(AT, &alpha, x, &beta, y); 
  return 0;
}

static int psolve(double *x, double *b)
{ 
  psolvemx(A,&pivot,LU,D,x,b);
  return 0; 
}

static int psolvetrans(double *x, double *b)
{ 
  psolvemx(AT,&pivotT,LUT,D,x,b);  
  return 0;
}

int estiva_bicgsolver(void* Ap, double* x, double* b)
{
  long i, n, ldw, iter, info;
  static double *work, resid;

  n     = dim1(b);
  ldw   = n;
  iter  = 100 * n;
  resid = 0.0000001;

  ary1(D,n+1);
  ary1(work,n*6+1);

  slimupmx(A,Ap);
  transmx(AT,A);

  if ( !strcmp(getop("-precond"),"ILU") ) {
    clonemx(LU,A);   
    ILU(&pivot,LU);
    clonemx(LUT,AT); 
    ILU(&pivotT,LUT);
  }
  
  forall(0,i,n-1)
    if (mx(A,i+1,i+1) == 0.0) {
      D[i] = 1.0;
    }
    else {
      D[i] = 1.0/mx(A,i+1,i+1);
    }
  
  forall(0,i,n-1) x[i] = b[i];
  
  bicg_(&n, b, x, &work[1],
	 &ldw, &iter, &resid, matvec, matvectrans, psolve, psolvetrans, &info);
  
  if ( defop("-v") ) printf("iter = %ld\n",iter);

  return iter;
}


static int addvec(double da, double *dx,  double *dy)
{
  long i;
  forall(0,i,dim1(dy)) dy[i] += da * dx[i];
  return 0;
}

static double dotvec(double *dx, double *dy)
{
  double tmp = 0.0;
  long i;
  forall(0,i,dim1(dy)) tmp += dy[i] * dx[i];  
  return tmp;
}

static double L2(double *dx)
{
  double sum = 0.0;
  long i;
  forall(0,i,dim1(dx)) sum += dx[i]*dx[i];
  return sqrt(sum);
}

static void cpy(double *src, double *dst)
{
  long i;
  forall(0,i,dim1(dst)) dst[i] = src[i];
}

static int bicg_(long *n, double *bp, double *xp, double *work, 
		 long *ldw, long *iter, double *resid, 
		 int (*matvec)(), int (*matvectrans)(), 
		 int (*psolve)(), int (*psolvetrans)(), long *info)
{
  static double *p, *ptld, *q, *qtld, *r, *rtld, *z, *ztld, *x, *b;
  static double alpha, beta, bnrm2, rho, rho1, rhotol, tol;
  static long   maxit;
  long   work_dim1,i;
  
  work_dim1   = *ldw;
  *info       = 0;
  maxit       = *iter;
  tol         = *resid;
  *iter       = 0;

  ary1(r,work_dim1);
  ary1(rtld,work_dim1);
  ary1(z,work_dim1);
  ary1(ztld,work_dim1);
  ary1(p,work_dim1);
  ary1(ptld,work_dim1);
  ary1(q,work_dim1);    
  ary1(qtld,work_dim1);

  ary1(x,work_dim1+1);
  ary1(b,work_dim1+1);

  forall(0,i,work_dim1) x[i] = xp[i+1];
  forall(0,i,work_dim1) b[i] = bp[i+1];
  
  rhotol = 0.00000000000000000000000000000001232595164407830945955825883254353483864385054857848444953560829162597656250;
  
  cpy(b,r);
  if (L2(x) != 0.0) {
    matvec(-1.0, x, 1.0, r);
    if (L2(r) <= tol) {
      return 0;
    }
  }
  cpy(r,rtld);
  bnrm2 = L2(b);
  if (bnrm2 == 0.0) bnrm2 = 1.0;
  while ( (*iter)++ < maxit ) {
    psolve(z, r);
    psolvetrans(ztld, rtld);
    rho = dotvec(z, rtld);
    if ( fabsl(rho) < rhotol ) {
      *info = -10;
      return 0;
    }
    if (*iter > 1) {
      beta = rho / rho1;
      addvec(beta, p, z);
      addvec(beta, ptld, ztld);
      cpy(z,p);
      cpy(ztld,ptld);
    } else {
      cpy(z,p);
      cpy(ztld,ptld);
    }
    matvec(1.0, p, 0.0, q);
    matvectrans(1.0, ptld, 0.0, qtld);
    alpha = rho / dotvec(ptld,q);
    addvec( alpha, p, x);
    addvec(-alpha, q, r);
    *resid = L2(r) / bnrm2;
    if (*resid <= tol) break;
    addvec(-alpha, qtld,rtld);
    rho1 = rho;
  }

  forall(0,i,work_dim1) xp[i] = x[i-1];
  forall(0,i,work_dim1) bp[i] = b[i-1];
  
  *info = 1;
  return 0;
}
