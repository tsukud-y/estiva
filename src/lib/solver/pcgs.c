#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"
#include "estiva/eblas.h"

int estiva_pcgssolver(void *Apointer, double *xk, double *b)
{
  static double *ek, *ek1, *hk1, *pk, *pk1, *q, *r0, *rk, *rk1, *w, *xk1, *xk1_xk;
  MX *A;
  double c1, c2, c3, alphak, betak;
  long   k, n;

  A = Apointer;
  n = A->n;
  setAmx(A);
  ILUdecomp(A);

  setveclength(n+1);
  if ( defop("-adjust") ) {
    setveclength(n);
    b=&b[1];
    xk=&xk[1];
  }
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

  matvecvec(A,1.0,xk,0.0,q);
  addformula( rk, '=', b, '-', 1.0,q);
  psolvevec(A,rk);
  cpvec(rk,r0);
  cpvec(rk,ek);
  cpvec(rk,pk);
  c1 = dotvec(r0,r0);

  forall (1, k, n) {
    phase(1); {
      matvecvec(A, 1.0, pk, 0.0, q);
      psolvevec(A,q);
      c2 = dotvec(q,r0);
      alphak = c1 / c2;
    }
    phase(2); {
      addformula( hk1, '=', ek,'-',alphak,q);
      addformula( w,   '=', ek,'+',1.0,hk1); 
    }
    phase(3); {
      matvecvec(A, 1.0, w, 0.0, q);
      psolvevec(A,q);
    }
    phase(4); {
      addformula( rk1, '=', rk,'-',alphak,q);
      addformula( xk1, '=', xk,'+',alphak,w);
      c3 = dotvec(rk1,r0);
    }
    phase(5); {
      addformula( xk1_xk, '=', xk1,'-',1.0,xk);
      if ( L2(xk1_xk)/L2(xk) < epsilon() && psc98condition(xk1,b) ) {
	  cpvec(xk1,xk);
	  return success(k);
      }
    }
    phase(6); {
      betak = c3 / c1;
      c1 = c3;
      addformula( ek1, '=', rk1,'+',betak,hk1);
      addformula( pk1, '=', ek1,'+',betak,addformula( w, '=', hk1,'+',betak,pk));
    }
    cpvec(ek1,ek);
    cpvec(pk1,pk);
    cpvec(rk1,rk);
    cpvec(xk1,xk);
  }
  return 1;
}
