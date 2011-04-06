#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"

int estiva_bicgsolver(void *Apointer, double *x, double *b)
{
  static MX *A, *AT;
  static double *p, *ptld, *q, *qtld, *r, *rtld, *z, *ztld, *dd;
  double alpha, beta, bnrm2, rho, rho1=1.0;
  long   n, itr;

  A = Apointer;

  transmx(AT,A);
  ILUdecomp(A);
  n = A->n;
  setveclength(n+1);
  if( defop("-adjust") ) {
    setveclength(n);
    x = &x[1];
    b = &b[1];
  }
  ary1(r   ,n+1);
  ary1(rtld,n+1);
  ary1(z   ,n+1);
  ary1(ztld,n+1);
  ary1(p   ,n+1);
  ary1(ptld,n+1);
  ary1(q   ,n+1);    
  ary1(qtld,n+1);
  ary1(dd,  n+1);

  cpvec(b,x);
  
  cpvec(b,r);
  if (L2(x) != 0.0) {
    matvecvec(A,-1.0, x, 1.0, r);
    if (L2(r) <= epsilon()) {
      return 0;
    }
  }
  cpvec(r,rtld);
  bnrm2 = L2(b);
  if (bnrm2 == 0.0) bnrm2 = 1.0;

  forall (1, itr, n) {
    cpvec(r,z);
    psolvevec(A,z);
    cpvec(rtld,ztld);
    psolvevec(AT,ztld);
    rho = dotvec(z, rtld);
    if ( fabs(rho) < 1.2e-31 ) break;
    if (itr == 1) {
      cpvec(z,p);
      cpvec(ztld,ptld);
    } else {
      beta = rho / rho1;
      addformula(z, '=', z, '+', beta,p );
      cpvec(z,p);
      addformula(ztld, '=', ztld, '+', beta,ptld );
      cpvec(ztld,ptld);
    }
    matvecvec(A, 1.0, p, 0.0, q);
    matvecvec(AT,1.0, ptld, 0.0, qtld);
    alpha = rho / dotvec(ptld,q);
    addformula(    x, '=', x,    '+', alpha,p    );
    addformula(    r, '=', r,    '-', alpha,q    );
    addformula( rtld, '=', rtld, '-', alpha,qtld );

    if (L2(r)/bnrm2 <= epsilon() && stopcondition(A,x,b)) 
      return success(itr);

    rho1 = rho;
  }
  return 1;
}
