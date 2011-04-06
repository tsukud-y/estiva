#include <math.h>
#include <stdio.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"

int estiva_cgssolver(void *Apointer, double *x, double *b)
{
  static double *dd, *p, *phat, *q, *qhat, *r, *rtld, *tmp, *u, *uhat, *vhat;
  MX *A;
  double alpha, beta, bnrm2, rho, rho1=1.0;
  long   itr, n;

  A = Apointer;
  n = A->n;
  ILUdecomp(A);

  setveclength(n+1);
  if ( defop("-adjust") ) {
    setveclength(n);
    b=&b[1];
    x=&x[1];
  }
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

  cpvec(b,x);
  cpvec(b, r);
  if ( L2(x) != 0.0 ) {
    matvecvec(A,-1.0, x, 1.0, r);
    if (L2(r) <= 1.0e-7) return 0;
  }
  bnrm2 = L2(b);
  if (bnrm2 == 0.0) bnrm2 = 1.0;
  cpvec(r, rtld);

  for (itr = 1; itr < n;  itr++) {
    rho = dotvec(rtld,r);
    if (fabs(rho) < 1.2e-31) break;
    if ( itr == 1 ) {
      cpvec(r,u);
      cpvec(u,p);
    } else {
      beta = rho / rho1;
      addformula( u, '=', r, '+', beta, q);
      addformula( p, '=', u, '+', beta, addformula( tmp, '=',q,'+',beta,p));
    }
    cpvec(p,phat);
    psolvevec(A,phat);
    matvecvec(A,1.0, phat, 0.0, vhat );

    alpha = rho / dotvec(rtld, vhat);
    addformula( q, '=', u, '-', alpha, vhat);
    cpvec(addformula(phat, '=', u, '+', 1.0, q),uhat);
    psolvevec(A,uhat);
    
    addformula(x, '=', x, '+', alpha, uhat);
    
    matvecvec(A,1.0, uhat, 0.0, qhat );
    addformula(r, '=', r, '-',alpha,qhat);
    
    if ( L2(r) / bnrm2 <= epsilon() && stopcondition(A,x,b) ) 
      return success(itr);
    rho1 = rho;
  }
  return 1;
}
