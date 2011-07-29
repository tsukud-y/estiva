#pragma OPENCL EXTENSION all: enable


static long get_vec_length(long n)
{
  __local long vec_length;
  if ( n == 0 ) return vec_length;
  vec_length = n;
  return  0;
}

static __global double *get_z_vec(__global double *z)
{
  __local long n;
  if ( z == 0 ) return (__global double *)n;
  n = (long)z;
  return 0;
}

static __global double *get_t_vec(__global double *z)
{
  __local long n;
  if ( z == 0 ) return (__global double *)n;
  n = (long)z;
  return 0;
}

static double dotgpu(__global double *x, __global double *y)
{
  __global double *z;
  long id, np, n, start, end, i;

  z     = get_z_vec(0);
  id    = get_global_id(0);	
  np    = get_global_size(0);
  n     = get_vec_length(0);
  start =  1+ id*n/np;
  end   = (1+id)*n/np;

  for (z[id] = 0.0, i=start; i<=end; i++) z[id] += x[i-1]*y[i-1];
  barrier(CLK_GLOBAL_MEM_FENCE);

  if(id == 0) for ( i=1; i<=np; i++) z[0] += z[i];
  barrier(CLK_GLOBAL_MEM_FENCE);

  return z[0];
}

__kernel void dotParallel(__global double *x, __global double *y, __global double *z, long n)
{
  n++;
  get_vec_length(n);
  get_z_vec(z);
  
  z[0] =  dotgpu(x,y);
}

static void cpgpu(__global double *x, __global double *y)
{
  long id, np, n, start, end, i;

  id    = get_global_id(0);	
  np    = get_global_size(0);
  n     = get_vec_length(0);
  start =  1+ id*n/np;
  end   = (1+id)*n/np;
  
  for (i=start; i<=end; i++) y[i-1] = x[i-1];
}

__kernel void cpParallel(
	 __global double* x,
	 __global double* y,
	 long n)
{
  n++;
  get_vec_length(n);
  cpgpu(x,y);
}


static void addgpu(double da, __global double *x, __global double *y)
{
  long id, np, n, start, end, i;

  id    = get_global_id(0);	
  np    = get_global_size(0);
  n     = get_vec_length(0);
  start =  1+ id*n/np;
  end   = (1+id)*n/np;
  
  for (i=start; i<=end; i++) y[i-1] += da * x[i-1];
}

__kernel void addParallel(double da, __global double* x, __global double* y, long n)
{
  n++;
  get_vec_length(n);
  addgpu(da,x,y);
}

static double L2gpu(__global double *x)
{
  return sqrt(dotgpu(x,x));
}

__kernel void L2Parallel(__global double *x, __global double *z, long n)
{
  n++;
  
  get_vec_length(n);
  get_z_vec(z);
  
  z[0] = L2gpu(x);
}

typedef struct vcr_s {
  __global double *val;
  __global long *col_ind;
  __global long *row_ptr;
  int set;
} vcr_struct;

static vcr_struct get_vcr_struct(vcr_struct a)
{
  __local vcr_struct A;

  if ( a.set == 1) {
    A.val     = a.val;
    A.col_ind = a.col_ind;
    A.row_ptr = a.row_ptr;
    A.set     = 0;
  }
  return A;
}

static void mulgpu(vcr_struct A, __global double *x, __global double *y)
{
  long i, j, start, goal, n;

  n     = get_vec_length(0);
  start = 1+(get_global_id(0))*n/get_global_size(0);
  goal  = (get_global_id(0)+1)*n/get_global_size(0);

  barrier(CLK_GLOBAL_MEM_FENCE);
  for (i=start; i<=goal; i++) {
    y[i-1] = 0;
    for (j=A.row_ptr[i]; j<A.row_ptr[i+1]; j++)
      y[i-1] += A.val[j]*x[A.col_ind[j]-1];
  } 
}


static void matvecgpu(vcr_struct A, double alpha, __global double *x, double beta, __global double *y)
{
  __global double *t;

  t = get_t_vec(0);
  mulgpu(A,x,t);
  addgpu(alpha-1.0,t,t);
  addgpu(beta,y,t);
  cpgpu(t,y);
}

__kernel void dataParallel(
			   __global double* val, 
			   __global long* col_ind, 
			   __global long* row_ptr,
			   double alpha,
			   __global double* x,
			   double beta,
			   __global double* y,
			   __global double* t,
			   long n)
{ 
  vcr_struct A;

  A.val     = val;
  A.col_ind = col_ind;
  A.row_ptr = row_ptr;
  A.set     = 1;
  get_vcr_struct(A);
  A = get_vcr_struct(A);

  get_vec_length(n);
  get_t_vec(t);

  matvecgpu(A,alpha,x,beta,y);
}


__kernel void cg_iterParallel( __global double *val, __global long *col_ind, __global long *row_ptr,
			       __global double *x,
			       __global double *pvec,
			       __global double *qvec,
			       __global double *rvec,
			       __global double *tvec,
			       __global double *zvec,
			       double bnrm2, double rho1, double tol, long maxit,
			       long n
			      )
{
  vcr_struct A;
  double alpha, beta, rho;
  long *iter, iterr = 0;
  double *resid, residr;
  
  iter = &iterr;
  resid = &residr;

  A.val     = val;
  A.col_ind = col_ind;
  A.row_ptr = row_ptr;
  A.set     = 1;
  get_vcr_struct(A);
  A = get_vcr_struct(A);

  get_vec_length(n);
  get_t_vec(tvec);
  get_z_vec(tvec);

 L10:
  /*        Perform Preconditioned Conjugate Gradient iteration. */
  ++(*iter);
  /*        Preconditioner Solve. */

  //(*psolve)(zvec, rvec);                                                                                                     
  cpgpu(rvec,zvec);
  
  rho = dotgpu(zvec,rvec);
  /*        Compute direction vector P. */
  if (*iter > 1) {
    beta = rho / rho1;
    addgpu(beta, pvec, zvec);
    cpgpu(zvec, pvec);
  } else {
    cpgpu(zvec, pvec);
  }
  /*        Compute scalar ALPHA (save A*P to Q). */

  //matvecgpu(A, 1.0, pvec, 0.0, qvec);
  mulgpu(A,pvec,qvec);
  alpha = rho / dotgpu( pvec, qvec);

  /*        Compute current solution vector X. */
  addgpu(alpha, pvec,  x);

  /*        Compute residual vector R, find norm, */
  /*        then check for tolerance. */
  addgpu(-alpha, qvec, rvec);
  *resid = L2gpu(rvec) / bnrm2;

  //if (*resid <= tol && psc98condition(x,b)) goto L30;                                                                        
  if (*resid <= tol && (*iter) >= 343) goto L30;
  if (*iter == maxit) goto L20;

  rho1 = rho;

  goto L10;

 L20:
  /*     Iteration fails. */
  return ;

 L30:
  /*     Iteration successful; return. */
  return ;
  /*     End of CG */
}
