#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/std.h"
#include "estiva/vec.h"
#include "estiva/op.h"
#include "estiva/eblas.h"
#include "estiva/solver.h"
#include "estiva/matvecmpi2.h"
#include <CL/cl.h>
#include "tocuda.h"
#include "Rf5.h"

static int cg_();

int estiva_cgsolver2(void *A, double *x, double *b)
{
  static double *b1, *x1;
  static double *work;
  long n, ldw, iter, info, i;
  double resid = 1.0e-7;


  mxtocuda(A);


  if ( !symcheckmx(A) ) {
    fprintf(stderr,"matrix is not symmetric\n");
    return 1;
  }


  setAmx(A);
  ldw = iter = n = dim1(b);

  setveclength(n);
  forall (0, i, n ) x[i] = b[i];
  ILUdecomp(A);

  ary1(b1,n);
  ary1(x1,n);

  forall (1,i,n) b1[i-1] = b[i];
  forall (1,i,n) x1[i-1] = x[i];
  

  cg_(A,&n,b1,x1,work,&ldw,&iter,&resid,estiva_matvec,estiva_psolve,&info);
  if (defop("-v")) fprintf(stderr, "cg iter = %ld\n",iter);
  
  forall (1,i,n) x[i] = x1[i-1];




  return 0;
}


/* CG.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* Table of constant values */

static long c__1 = 1;
static double c_b5 = -1.;
static double c_b6 = 1.;
static double c_b20 = 0.;

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
*  CG solves the linear system Ax = b using the
*  Conjugate Gradient iterative method with preconditioning.
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
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,4).
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
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2
*  ============================================================ */
#include "Rf4.h"

void estiva_matvecParallel(MX *A, double alpha, double *x, double beta, double *y)
{
  static double *t;
  CRSMATRIX *a;
  long i, n = A->n; 
  size_t global_item_size = 1;
  size_t local_item_size = 1;
  cl_kernel kernel = estiva_clmulCRSmatrixvec;
  cl_int ret;
  ary1(t,n);

  a = *(CRSMATRIX**)estiva_std_f4(A);

  if (defop("-gpu")) global_item_size= atoi(getop("-gpu"));
  local_item_size = global_item_size;

  ary1tocuda(t);
  ary1tocuda(x);
  ary1tocuda(y);

  ret = clSetKernelArg(kernel, 0, ary1arg(a->val));
  ret = clSetKernelArg(kernel, 1, ary1arg(a->col_ind));
  ret = clSetKernelArg(kernel, 2, ary1arg(a->row_ptr));
  ret = clSetKernelArg(kernel, 3, ary0arg(alpha));
  ret = clSetKernelArg(kernel, 4, ary1arg(x));
  ret = clSetKernelArg(kernel, 5, ary0arg(beta));
  ret = clSetKernelArg(kernel, 6, ary1arg(y));
  ret = clSetKernelArg(kernel, 7, ary1arg(t));
  ret = clSetKernelArg(kernel, 8, ary0arg(n));

  ret = clEnqueueNDRangeKernel(clcommand_queue(), kernel, 1, NULL,
                               &global_item_size, &local_item_size, 0, NULL, NULL);
  if (ret) abort();
  ary1fromcuda(y);
}

extern double estiva_L2Parallel(double *x);
extern double estiva_dotParallel(double *x, double  *y);

static int cg_iter(MX *A, double *x, double *pvec, double *qvec, double *rvec, double *zvec,
		   double bnrm2, double rho1, double tol, long maxit, double *resid, long *iter, double *b)
{
  double alpha, beta, rho;

  static double *tvec;
  CRSMATRIX *a;
  long  n = A->n;
  size_t global_item_size = 1;
  size_t local_item_size = 1;
  cl_kernel kernel = estiva_clcg_iter;
  cl_int ret = 0;

  ary1(tvec,n);
  
  a = *(CRSMATRIX**)estiva_std_f4(A);

  if (defop("-gpu")) global_item_size= atoi(getop("-gpu"));
  local_item_size = global_item_size;

  ary1tocuda(x);
  ary1tocuda(pvec);
  ary1tocuda(qvec);
  ary1tocuda(rvec);
  ary1tocuda(tvec);
  ary1tocuda(zvec);
  ary1tocuda(resid);
  ary1tocuda(iter);

  clSetKernelArg(kernel, 0, ary1arg(a->val));
  clSetKernelArg(kernel, 1, ary1arg(a->col_ind));
  clSetKernelArg(kernel, 2, ary1arg(a->row_ptr));
  clSetKernelArg(kernel, 3, ary1arg(x));
  clSetKernelArg(kernel, 4, ary1arg(pvec));
  clSetKernelArg(kernel, 5, ary1arg(qvec));
  clSetKernelArg(kernel, 6, ary1arg(rvec));
  clSetKernelArg(kernel, 7, ary1arg(tvec));
  clSetKernelArg(kernel, 8, ary1arg(zvec));
  clSetKernelArg(kernel, 9, ary0arg(bnrm2));
  clSetKernelArg(kernel,10, ary0arg(rho1));
  clSetKernelArg(kernel,11, ary0arg(tol));
  clSetKernelArg(kernel,12, ary0arg(maxit));
  clSetKernelArg(kernel,13, ary0arg(n));

  if ( !defop("-host") ) {
    ret = clEnqueueNDRangeKernel(clcommand_queue(), kernel, 1, NULL,
			   &global_item_size, &local_item_size, 0, NULL, NULL);
    ary1fromcuda(x);
    return 0;
  }
 L10:
  /*        Perform Preconditioned Conjugate Gradient iteration. */
  ++(*iter);
  printf("iter = %d\n",*iter);
  /*        Preconditioner Solve. */


  //(*psolve)(zvec, rvec);
  estiva_cpParallel(rvec,zvec);


  rho = estiva_dotParallel(zvec,rvec);

  /*        Compute direction vector P. */
  
  if (*iter > 1) {
    beta = rho / rho1;
    estiva_addParallel(beta, pvec, zvec);
    estiva_cpParallel(zvec, pvec);
  } else {
    estiva_cpParallel(zvec, pvec);
  }
  
  /*        Compute scalar ALPHA (save A*P to Q). */
  estiva_matvecParallel(A, 1.0, pvec, 0.0, qvec);


  
  alpha = rho / estiva_dotParallel( pvec, qvec);

  
  /*        Compute current solution vector X. */
  estiva_addParallel(alpha, pvec,  x);




  /*        Compute residual vector R, find norm, */
  /*        then check for tolerance. */
  estiva_addParallel(-alpha, qvec, rvec);

  *resid = estiva_L2Parallel(rvec) / bnrm2;
  

  if (*resid <= tol && psc98condition(x,b)) goto L30;
  //if (*resid <= tol && (*iter) >= 343) goto L30;
  if (*iter == maxit) goto L20;
  
  rho1 = rho;
  goto L10;
  
 L20:
  /*     Iteration fails. */
  return 1;
  
 L30:
  /*     Iteration successful; return. */
  return 0;
  /*     End of CG */
}


static int cg_(A, n, b, x, work, ldw, iter, resid, matvec, psolve, info)
     MX *A;
     long *n, *ldw, *iter, *info;
     double *b, *x, *work, *resid;
     int (*matvec) (), (*psolve) ();
{
    /* System generated locals */
    /* Local variables */
    static double *pvec, *qvec, *rvec, *zvec; 
    double beta, bnrm2, alpha, rho, tol, rho1;
    long  maxit;

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
    ary1(rvec,*n);
    ary1(zvec,*n);
    ary1(pvec,*n);
    ary1(qvec,*n);

    /*     Set initial residual. */
    estiva_cpParallel(b, rvec);
    if (estiva_L2Parallel(x) != 0.) {
      estiva_matvecParallel(A, c_b5, x, c_b6, rvec);
      if (estiva_L2Parallel(rvec) < tol)  goto L30;
    }

    bnrm2 = estiva_L2Parallel(b);

    if (bnrm2 == 0.) bnrm2 = 1.;
    *iter = 0;


    /*
    extern int cg_iter(long *iter, double *rvec, double *zvec, double rho, double *pvec, MX *A, double *qvec,
			 double *x, double bnrm2, double *regid, double tol, long maxit);
    */


    return cg_iter(A, x, pvec, qvec, rvec, zvec,
		   bnrm2, rho1, tol, maxit, resid, iter,b);
 L10:
    /*        Perform Preconditioned Conjugate Gradient iteration. */
    ++(*iter);
    printf("iter = %d\n", *iter);
    /*        Preconditioner Solve. */

    //(*psolve)(zvec, rvec);
    estiva_cpParallel(rvec,zvec);

    rho = estiva_dotParallel(zvec,rvec);

    /*        Compute direction vector P. */

    if (*iter > 1) {
	beta = rho / rho1;
	estiva_addParallel(beta, pvec, zvec);
	estiva_cpParallel(zvec, pvec);
    } else {
	estiva_cpParallel(zvec, pvec);
    }

    /*        Compute scalar ALPHA (save A*P to Q). */
    estiva_matvecParallel(A, c_b6, pvec, c_b20, qvec);
    alpha = rho / estiva_dotParallel( pvec, qvec);

    /*        Compute current solution vector X. */
    estiva_addParallel(alpha, pvec,  x);

    /*        Compute residual vector R, find norm, */
    /*        then check for tolerance. */

    estiva_addParallel(-alpha, qvec, rvec);

    *resid = estiva_L2Parallel(rvec) / bnrm2;


    //if (*resid <= tol && psc98condition(x,b)) goto L30;
    if (*resid <= tol && (*iter) >= 343) goto L30;
    if (*iter == maxit) goto L20;

    rho1 = rho;
    goto L10;

 L20:
    /*     Iteration fails. */
    *info = 1;
    return 0;
    
 L30:
    /*     Iteration successful; return. */
    return 0;
    /*     End of CG */
}
