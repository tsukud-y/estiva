#include <stdio.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <math.h>
#include <estiva/precond.h>
#include <estiva/solver.h>
#include <estiva/op.h>
#include "estiva/std.h"


static MX *A, *AT, *LU, *LUT;
static double *D;
static CRS pivot, pivotT;
static int cgs_();


static int matvec(double alpha, double *x, double beta, double *y)
{ 
  matvecmx(A, &alpha, x, &beta, y); 
  return 0;
}

static int psolve(double *x, double *b)
{
  psolvemx(A,&pivot,LU,D,x,b); 
  return 0; 
}


int estiva_cgssolver(void* pA, double* x, double* b)
{
   long int n,ldw,iter,info, i; 
   static double *work, resid;

   A = pA;
   transmx(AT,A);


   if ( !strcmp(getop("-precond"),"ILU") ) {
     LU = NULL;
     clonemx(LU,A);
     ILU((void*)&pivot,LU);

     LUT = NULL;
     clonemx(LUT,AT);
     ILU((void*)&pivotT,LUT);
   }

   n = dim1(b);

   ary1(D,n+1);
   for(i=0;i<n;i++) 
     if (mx(A,i+1,i+1) == 0.0) D[i] = 1.0;
     else                      D[i] = 1.0/mx(A,i+1,i+1);

   ary1(work,n*7+1);
   ldw = n;
   iter = 100 * n;
   resid = 0.0000001;

   for ( i=0; i<n; i++ ) x[i] = b[i];

   cgs_(&n, &b[1], &x[1], &work[1],
	 &ldw, &iter, &resid, matvec, psolve, &info);



   printf("iter = %ld\n",iter);

   return iter;
}

/********************************* CGS.c ************************************/
/* CGS.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  Purpose
*  =======
*
*  CGS solves the linear system Ax = b using the
*  Conjugate Gradient Squared iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*         The preconditioner is passed into the routine in a common block
*
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If RHO become smaller than some tolerance,
*                   the program will terminate. Here we check
*                   against tolerance BREAKTOL.
*
*                   -10: RHO < BREAKTOL: RHO and RTLD have become
*                                        orthogonal.
*
*  BLAS CALLS:    DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ============================================================ */


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

static void scal(double da, double *dx)
{
  long i;
  forall(0,i,dim1(dx)) dx[i] *= da;
}

static int cgs_(n, b, x, work, ldw, iter, resid, matvec, psolve, info)
     long *n, *ldw, *iter, *info;
     double *b, *x, *work, *resid;
     int (*matvec) (), (*psolve) ();
{

  double bnrm2;
  static double *rtld;
  static double *phat, *qhat, *uhat, *vhat;
  static double *p, *q, *r, *u;
  static double rhotol, rho, tol, rho1;
  static double *x0, *b0;
  double alpha, beta;
  long i, maxit, work_dim1;  

  /* Parameter adjustments */
  work_dim1 = *ldw;
  /* Executable Statements */
  *info = 0;
  /* Test the input parameters. */
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
  /*  Alias workspace columns. */
  ary1(r    ,work_dim1);
  ary1(rtld ,work_dim1);
  ary1(p    ,work_dim1);
  ary1(phat ,work_dim1);
  ary1(q    ,work_dim1);
  ary1(qhat ,work_dim1);
  ary1(u    ,work_dim1);
  ary1(uhat ,work_dim1);
  ary1(vhat ,work_dim1);
  ary1(x0,   work_dim1);
  ary1(b0,   work_dim1);
  forall(0,i,dim1(x0)) x0[i] = x[i];
  forall(0,i,dim1(b0)) b0[i] = b[i];
  /* Set breakdown tolerance parameter. */
  rhotol =  0.00000000000000000000000000000001232595164407830945955825883254353483864385054857848444953560829162597656250;
  /* Set initial residual. */
  cpy(b0, r);
  if (L2(x0) != 0.0) {
    matvec(-1.0, x0, 1.0, r);
    if (L2(r) <= tol) {
      return 0;
    }
  }
  bnrm2 = L2(b0);
  if (bnrm2 == 0.0) {
    bnrm2 = 1.0;
  }
  /* Choose RTLD such that initially, (R,RTLD) = RHO is not equal to 0. */
  /* Here we choose RTLD = R. */
  cpy(r, rtld);

  for (*iter = 1; ; (*iter)++) {
    /* Perform Conjugate Gradient Squared iteration. */
    rho = dotvec(rtld,r);
    if (fabsl(rho) < rhotol) {
      /* Set breakdown flag. */
      if (fabsl(rho) < rhotol) {
	*info = -10;
      }
      break;
    }
    /* Compute direction vectors U and P. */
    if (*iter > 1) {
      /* Compute U. */
      beta = rho / rho1;
      cpy(r,u);
      addvec(beta,q,u);
      /* Compute P. */
      /* Computing 2nd power */
      scal(beta*beta, p);
      addvec(beta, q, p);
      addvec(1.0, u, p);
    } else {
      cpy(r,u);
      cpy(u,p);
    }
    /* Compute direction adjusting scalar ALPHA. */
    psolve(phat, p);
    matvec(1.0, phat, 0.0, vhat );
    alpha = rho / dotvec(rtld, vhat);
    cpy(u, q);
    addvec(-alpha, vhat, q);
    /* Compute direction adjusting vector UHAT. */
    /* PHAT is being used as temporary storage here. */
    cpy(q, phat);
    addvec(1.0, u, phat);
    psolve(uhat, phat);
    /* Compute new solution approximation vector X. */
    addvec(alpha, uhat, x0);
    /* Compute residual R and check for tolerance. */
    matvec(1.0, uhat, 0.0, qhat );
    addvec(-alpha, qhat,  r);
    *resid = L2(r) / bnrm2;
    if (*resid <= tol) {
      /* Iteration successful; return. */
      break;
    }
    if (*iter == maxit) {
      *info = 1;
      /* Iteration fails. */
      break;
    }
    rho1 = rho;
  }
  forall(0,i,dim1(x0)) x[i] = x0[i];
  forall(0,i,dim1(b0)) b[i] = b0[i];
  return 0;
  /*  End of CGS */
}

