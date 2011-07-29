#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
#include <estiva/ary.h>
#include <estiva/std.h>
#include <estiva/mx.h>
#include <estiva/op.h>
#include "tocuda.h"
#include "Rf4.h"
#include "Rf5.h"
#include "estiva/vec.h"

#define MAX_SOURCE_SIZE (0x1000000)

char            *estiva_tocuda_source_str    = NULL;
cl_command_queue estiva_tocuda_command_queue = NULL;
cl_context       estiva_tocuda_context       = NULL;
cl_program       estiva_tocuda_program       = NULL;
cl_device_id     estiva_tocuda_device_id     = NULL;
cl_kernel        estiva_clmulCRSmatrixvec = NULL;
cl_kernel        estiva_cldotvec = NULL;
cl_kernel        estiva_clL2Parallel = NULL;
cl_kernel        estiva_clcp = NULL;
cl_kernel        estiva_cladd = NULL;
cl_kernel        estiva_clcg_iter = NULL;

void  estiva_tocuda(const char *fileName)
{
  static int tocuda_init = 0;
  cl_platform_id platform_id = NULL;
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  size_t source_size;
  FILE *fp;
  cl_int ret;

  if (tocuda_init != 0) return;
  tocuda_init = 1;
  
  /* カーネルを含むソースコードをロード */
  fp = fopen(fileName, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }
  clsource_str() = (char *)malloc(MAX_SOURCE_SIZE);
  source_size = fread( clsource_str(), 1, MAX_SOURCE_SIZE, fp );
  fclose( fp );

  /* プラットフォーム・デバイスの情報の取得 */
  ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &cldevice_id(), &ret_num_devices);
      
  /* OpenCLコンテキストの作成 */
  clcontext() = clCreateContext( NULL, 1, &cldevice_id(), NULL, NULL, &ret);
  
  /* コマンドキューの作成 */
  clcommand_queue() = clCreateCommandQueue(clcontext(), cldevice_id(), 0, &ret);
  
  /* 読み込んだソースからカーネルプログラムを作成 */
  clprogram() = clCreateProgramWithSource(clcontext(), 1, (const char **)&clsource_str(), (const size_t *)&source_size, &ret);
  ret     = clBuildProgram(clprogram(), 1, &cldevice_id(), NULL, NULL, NULL);
  free(clsource_str());

  if (ret != CL_SUCCESS) {
    char buffer[1024];
    clGetProgramBuildInfo(clprogram(), cldevice_id(), CL_PROGRAM_BUILD_LOG, 1024, buffer, NULL);
    puts(buffer);
    abort();
  }
}


void estiva_tocuda_Finalize(void)
{
  cl_int ret;

  ret = clReleaseKernel(estiva_clmulCRSmatrixvec);
  if (ret) abort();
  ret = clFlush(clcommand_queue());
  if (ret) abort();
  ret = clFinish(clcommand_queue());
  if (ret) abort();
  ret = clReleaseProgram(clprogram());
  if (ret) abort();
  ret = clReleaseCommandQueue(clcommand_queue());
  if (ret) abort();
  ret = clReleaseContext(clcontext());
  if (ret) abort();
}

cl_int estiva_ary1tocuda(void *a)
{
  cl_int ret;
  if ( estiva_std_f5(a) == NULL ) {R5(a,calloc(1,sizeof(cl_mem)));
  *(cl_mem*)estiva_std_f5(a) = clCreateBuffer(clcontext(), CL_MEM_READ_WRITE, dimdim(a), NULL, &ret);
  }

  ret = clEnqueueWriteBuffer(clcommand_queue(),  (*(cl_mem*)estiva_std_f5(a)), CL_TRUE, 0, dimdim(a), a, 0, NULL, NULL );


  return ret;
}



void estiva_initCRSmatrix(CRSMATRIX **Ap)
{
  CRSMATRIX *T;
  if (*Ap == NULL) ary1(*Ap,1);
  T = *Ap;
  T->val = NULL;
  T->col_ind = NULL;
  T->row_ptr = NULL;
}


void estiva_mx2CRSmatrix(MX *A, CRSMATRIX *B)
{
  long i, j, n, N=0;

  fornonzeromx(A,i,j) N++;
  n = A->n;

  ary1(B->val,N+1);
  ary1(B->col_ind,N+1);
  ary1(B->row_ptr,n+2);

  N = 1;
  fornonzeromx(A,i,j) {
    B->val[N] = mx(A,i,j);
    B->col_ind[N] = j;
    N++;
  }


  //forall(1,i,dim1(B->val)) printf("%.0lf ",B->val[i]);                                               
  //printf("\n");                                                                                      

  //forall(1,i,dim1(B->col_ind)) printf("%ld ",B->col_ind[i]);                                         
  //printf("\n");                                                                                      

  N = 0;
  fornonzeromx(A,i,j) {
    if ( i == N ) B->row_ptr[i+1]++;
    else N++;
  }
  for (i=2; i<=dim1(B->row_ptr)+1; i++) {
    B->row_ptr[i-1]++;
  }
  for (i=2; i<=dim1(B->row_ptr); i++) {
    B->row_ptr[i] += B->row_ptr[i-1];
  }
  //forall(1,i,dim1(B->row_ptr)) printf("%ld\n",B->row_ptr[i]);                                        
}


void estiva_mulCRSmatrixvec(CRSMATRIX *A, double alpha, double *x, double beta, double *y, double *t)
/* OpenCLカーネル引数の設定 */
/* OpenCLカーネルをデータ並列で実行 */
{
  size_t global_item_size = 1;
  size_t local_item_size = 1;
  cl_kernel kernel = estiva_clmulCRSmatrixvec;
  cl_int ret;
  long n = dim1(A->row_ptr)-1;

  if (defop("-gpu")) global_item_size= atoi(getop("-gpu"));
  local_item_size = global_item_size;

  ary1tocuda(t);
  ary1tocuda(x);
  ary1tocuda(y);

  ret = clSetKernelArg(kernel, 0, ary1arg(A->val));
  ret = clSetKernelArg(kernel, 1, ary1arg(A->col_ind));
  ret = clSetKernelArg(kernel, 2, ary1arg(A->row_ptr));
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


double estiva_dotParallel(double *x, double *y)
{
  static double *z;
  cl_kernel kernel = estiva_cldotvec;
  size_t global_item_size = 1, local_item_size;
  long i, n = dim1(x);


  if (defop("-gpu")) global_item_size = atoi(getop("-gpu"));

  local_item_size = global_item_size;
  ary1(z,global_item_size+1);
  
  ary1tocuda(x);
  ary1tocuda(y);
  ary1tocuda(z);

  clSetKernelArg(kernel, 0, ary1arg(x));
  clSetKernelArg(kernel, 1, ary1arg(y));
  clSetKernelArg(kernel, 2, ary1arg(z));
  clSetKernelArg(kernel, 3, ary0arg(n));

  clEnqueueNDRangeKernel(clcommand_queue(), kernel, 1, NULL,  &global_item_size, &local_item_size, 0, NULL, NULL);

  ary1fromcuda(z);

  return z[0];
}

int estiva_cpParallel(double *x, double *y)
{
  cl_kernel kernel = estiva_clcp;
  size_t global_item_size = 1, local_item_size;
  long i, n = dim1(x);


  if (defop("-gpu")) global_item_size = atoi(getop("-gpu"));

  local_item_size = global_item_size;
  
  ary1tocuda(x);
  ary1tocuda(y);

  clSetKernelArg(kernel, 0, ary1arg(x));
  clSetKernelArg(kernel, 1, ary1arg(y));
  clSetKernelArg(kernel, 2, ary0arg(n));

  clEnqueueNDRangeKernel(clcommand_queue(), kernel, 1, NULL,  &global_item_size, &local_item_size, 0, NULL, NULL);

  ary1fromcuda(y);

  return 0;
}


int estiva_addParallel(double da, double *x, double *y)
{
  cl_kernel kernel = estiva_cladd;
  size_t global_item_size = 1, local_item_size;
  long i, n = dim1(x);


  if (defop("-gpu")) global_item_size = atoi(getop("-gpu"));

  local_item_size = global_item_size;
  
  ary1tocuda(x);
  ary1tocuda(y);

  clSetKernelArg(kernel, 0, ary0arg(da));
  clSetKernelArg(kernel, 1, ary1arg(x));
  clSetKernelArg(kernel, 2, ary1arg(y));
  clSetKernelArg(kernel, 3, ary0arg(n));

  clEnqueueNDRangeKernel(clcommand_queue(), kernel, 1, NULL,  &global_item_size, &local_item_size, 0, NULL, NULL);

  ary1fromcuda(y);

  return 0;
}


double estiva_L2Parallel(double *x)
{
  static double *z;
  cl_kernel kernel = estiva_clL2Parallel;
  size_t np = 1, lp = 1;
  long i, n = dim1(x);

  if (defop("-gpu")) np = atoi(getop("-gpu"));
  lp = np;
  ary1(z,np+1);
  ary1tocuda(x);
  ary1tocuda(z);

  clSetKernelArg(kernel, 0, ary1arg(x));
  clSetKernelArg(kernel, 1, ary1arg(z));
  clSetKernelArg(kernel, 2, ary0arg(n));

  clEnqueueNDRangeKernel(clcommand_queue(), kernel, 1, NULL,  &np, &lp, 0, NULL, NULL);

  ary1fromcuda(z);
  
  return z[0];
}



void estiva_mxtocuda_internal(MX *A, CRSMATRIX **Bp)
{
  cl_int ret;
  static CRSMATRIX *B;
  initCRSmatrix(*Bp);
  B = *Bp;
  mx2CRSmatrix(A,B);

  /* データ並列のOpenCLカーネルの作成 */
  tocuda("./dataParallel.cl");
  estiva_clmulCRSmatrixvec = clCreateKernel(clprogram(), "dataParallel", &ret);
  estiva_cldotvec = clCreateKernel(clprogram(),"dotParallel",&ret);
  estiva_clL2Parallel = clCreateKernel(clprogram(),"L2Parallel",&ret);
  estiva_clcp = clCreateKernel(clprogram(),"cpParallel",&ret);
  estiva_cladd = clCreateKernel(clprogram(),"addParallel",&ret);
  estiva_clcg_iter = clCreateKernel(clprogram(),"cg_iterParallel",&ret);

  /* バッファオブジェクトの作成 イニシャライズ */
  ret = ary1tocuda(B->val);
  ret = ary1tocuda(B->col_ind);
  ret = ary1tocuda(B->row_ptr);
  if (ret) abort();
}

void estiva_mxtocuda(MX *A)
{
  if ( estiva_std_f4(A) == NULL ) estiva_std_R4(A,calloc(1,sizeof(CRSMATRIX*)));
  estiva_mxtocuda_internal(A,estiva_std_f4(A));
}
