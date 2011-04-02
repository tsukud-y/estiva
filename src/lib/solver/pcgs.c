#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"

#if 0
int estiva_stopcondition(void *A, double *x, double *b);
int estiva_success(long k);
double estiva_epsilon(void);
int estiva_phase(int i);

#define stopcondition(A,x,b) estiva_stopcondition(A,x,b)
#define success(k)           estiva_success(k)
#define epsilon()            estiva_epsilon()
#define phase(i)             estiva_phase(i)

static double psc98(MX *A, double *x, double *b)
{
  static double *L;
  ary1(L,A->m);
  matvecvec(A,1.0,x,0.0,L);
  return Linf(addformula( L, '=', L, '-',1.0,b));
}

int estiva_stopcondition(void *Apointer, double *x, double *b)
{
  MX *A;
  A = Apointer;
  if ( defop("-psc98") ) {
    double norm;
    norm = psc98(A,x,b);
    if ( defop("-v") ) printf("%e\n",norm);
    if ( norm < 1.0e-8 ) return 1;
    else return 0;
  }	  
  return 1;
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
#endif

int estiva_pcgssolver(void *Apointer, double *xk, double *b)
{
  static double *ek, *ek1, *hk1, *pk, *pk1, *q, *r0, *rk, *rk1, *w, *xk1, *xk1_xk;
  MX *A;
  double c1, c2, c3, alphak, betak;
  long   k, n;

  A = Apointer;
  n = A->m;
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

  forall (1, k, max(n/10,1000)) {
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
      if ( L2(xk1_xk)/L2(xk) < epsilon() ) {
	if (stopcondition(A,xk1,b)) return success(k);
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
