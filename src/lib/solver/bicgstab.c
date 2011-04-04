#include <stdio.h>
#include <string.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"

static MX *A;
static int bicgstab_();

static int matvec(double *alpha, double *x, double *beta, double *y)
{
  matvecmx(A, alpha, x, beta, y); 
  return 0;
}

static int psolve(double *x, double *b)
{ 
  cpvec(b,x);
  psolvevec(A,x);
  return 0;
}

int estiva_bicgstabsolver(void *Apointer, double *x, double *b)
{
  static double *work;
  double resid;
  long int n, ldw, iter, info;
  
  A = Apointer;
  ILUdecomp(A);
  n = ldw = iter = A->m;
  setveclength(n);
  x = &x[1];
  b = &b[1];

  ary1(work,n*7);
  cpvec(b,x);
  
  resid = epsilon();
  
  bicgstab_(&n, b, x, work, &ldw, &iter, &resid, matvec, psolve, &info);
  
  return success(iter);
}

static double c_b5 = -1.;
static double c_b6 = 1.;
static double c_b25 = 0.;

static int bicgstab_(n, b, x, work, ldw, iter, resid, matvec, psolve, info)
   long *n, *ldw, *iter, *info;
   double *b, *x, *work, *resid;
   int (*matvec) (), (*psolve) ();
{
    /* System generated locals */
  long work_dim1;
    double d__1;

    /* Local variables */
    static double beta;

    static double omegatol, bnrm2;

    static double alpha;

    static double omega;

    static long maxit;

    static double rhotol, rho, tol, rho1;
    static double *r, *rtld, *p, *v, *t, *phat, *shat, *s;

    /* Parameter adjustments */
    work_dim1 = *ldw;


    --x;
    --b;

    /* Executable Statements */
    *info = 0;

/*     Test the input parameters. */

    if (*n < 0) {
	*info = -1;
    } else if (*ldw < max(1,*n)) {
	*info = -2;
    } else if (*iter <= 0) {
	*info = -3;
    }
    if (*info != 0) {
	return 0;
    }

    maxit = *iter;
    tol = *resid;

/*     Alias workspace columns. */

    ary1(r,work_dim1+1);
    ary1(rtld,work_dim1+1);
    ary1(p,work_dim1+1);
    ary1(v,work_dim1+1);
    ary1(t,work_dim1+1);
    ary1(phat,work_dim1+1);
    ary1(shat,work_dim1+1);
    ary1(s,work_dim1+1);

/*     Set parameter tolerances. */


    rhotol = omegatol = 1.2e-31;

/*     Set initial residual. */

    cpvec(b+1,r);
    if (L2(&x[1]) != 0.) {
	(*matvec)(&c_b5, &x[1], &c_b6, r);
	if (L2(r) <= tol && stopcondition(A,x,b)) {
	    goto L30;
	}
    }
    cpvec(r,rtld);

    bnrm2 = L2(&b[1]);
    if (bnrm2 == 0.) {
	bnrm2 = 1.;
    }

    *iter = 0;

L10:

/*     Perform BiConjugate Gradient Stabilized iteration. */

    ++(*iter);
    if ( defop("-v") ) fprintf(stderr,"iter = %d\n",*iter);

    rho = dotvec(rtld, r);
    if (fabsl(rho) < rhotol) {
	goto L25;
    }

/*        Compute vector P. */

    if (*iter > 1) {
	beta = rho / rho1 * (alpha / omega);
	d__1 = -omega;

	addvec(d__1, v, p);

	scalvec(beta, p);
	addvec(c_b6, r,p );
	
    } else {
      cpvec(r,p);
    }

/*        Compute direction adjusting vector PHAT and scalar ALPHA. */

    (*psolve)(phat   , p   );
    (*matvec)(&c_b6, phat   , &c_b25, v 
	    );
    alpha = rho / dotvec(rtld, v);

/*        Early check for tolerance. */

    d__1 = -alpha;
    //daxpy_(n, &d__1, v   , &c__1stab, r    , &c__1stab);
    addvec(d__1, v  , r);
    cpvec(r,s);
    if (L2(s) <= tol) {
      //daxpy_(n, &alpha, phat   , &c__1stab, &x[1], &c__1stab);
      addvec(alpha, phat, &x[1]);

	*resid = L2(s) / bnrm2;
	goto L30;
    } else {

/*           Compute stabilizer vector SHAT and scalar OMEGA. */

	(*psolve)(shat   , s   );
	(*matvec)(&c_b6, shat   , &c_b25, t 
		 );
	omega = dotvec(t, s) / dotvec(t, t);
/*           Compute new solution approximation vector X. */


	addvec(alpha,phat,&x[1]);
	addvec(omega, shat, &x[1]);

/*           Compute residual R, check for tolerance. */

	d__1 = -omega;

	addvec(d__1, t, r);

	*resid = L2(r) / bnrm2;
	if (*resid <= tol) {
	    goto L30;
	}
	if (*iter == maxit) {
	    goto L20;
	}
    }

    if (fabsl(omega) < omegatol) {
	goto L25;
    } else {
	rho1 = rho;
	goto L10;
    }

L20:

/*     Iteration fails. */

    *info = 1;
    return 0;

L25:

/*     Set breakdown flag. */

    if (fabsl(rho) < rhotol) {
	*info = -10;
    } else if (fabsl(omega) < omegatol) {
	*info = -11;
    }
    return 0;

L30:

/*     Iteration successful; return. */

    return 0;
/*     End of BICGSTAB */
}
