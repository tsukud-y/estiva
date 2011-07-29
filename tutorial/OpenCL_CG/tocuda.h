#ifndef _ESTIVA_TOCUDA_H_
#define _ESTIVA_TOCUDA_H_

typedef struct {
  double *val;
  long   *col_ind;
  long   *row_ptr;
} CRSMATRIX;

extern char            *estiva_tocuda_source_str;
extern cl_command_queue estiva_tocuda_command_queue;
extern cl_context       estiva_tocuda_context;
extern cl_program       estiva_tocuda_program; 
extern cl_device_id     estiva_tocuda_device_id;
extern void             estiva_tocuda(const char *fileName);
extern void             estiva_tocuda_Finalize(void);
extern cl_int estiva_ary1tocuda(void *a);
extern void estiva_initCRSmatrix(CRSMATRIX **Ap);
extern void estiva_mx2CRSmatrix(MX *A, CRSMATRIX *B);
extern cl_kernel estiva_clmulCRSmatrixvec;
extern cl_kernel estiva_clcg_iter;
extern void estiva_mulCRSmatrixvec(CRSMATRIX *A, double alpha, double *x, double beta, double *y, double *t);
extern void estiva_mxtocuda(MX *A);

#define clsource_str()    estiva_tocuda_source_str
#define clcommand_queue() estiva_tocuda_command_queue
#define clcontext()       estiva_tocuda_context
#define clprogram()       estiva_tocuda_program
#define cldevice_id()     estiva_tocuda_device_id
#define tocuda(fileName)  estiva_tocuda(fileName)
#define tocuda_Finalize() estiva_tocuda_Finalize()
#define ary1arg(a) sizeof(cl_mem), (void *)&ary1cl_mem(a)
#define ary0arg(n) sizeof(n),&(n)
#define dimdim(a) ((dim1(a)+1)*dim0(a))
#define ary1cl_mem(val) (*(cl_mem*)estiva_std_f5(val))
#define ary1tocuda(a) estiva_ary1tocuda(a)
#define clary1free(a) clReleaseMemObject(ary1cl_mem(a))
#define ary1fromcuda(y) clEnqueueReadBuffer(clcommand_queue(), ary1cl_mem(y), CL_TRUE, 0, dimdim(y), y, 0, NULL, NULL);
#define initCRSmatrix(A) estiva_initCRSmatrix(&(A))
#define mx2CRSmatrix(A,B) estiva_mx2CRSmatrix(A, B)
#define mulCRSmatrixvec(A,alpha,x,beta,y,t) estiva_mulCRSmatrixvec(A,alpha,x,beta,y,t)
#define mxtocuda(A) estiva_mxtocuda(A)

#endif
