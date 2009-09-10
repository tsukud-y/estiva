#ifndef _ARY_H_
#define _ARY_H_


#include <stdlib.h>
#define ary1(b,n_1) femlib_ary1((void *)&b,n_1,sizeof(*b))
#define ary2(A,m_1,n_1) femlib_ary2((void *)&A,m_1,n_1,sizeof(**A))

extern size_t dim0(void* v);
extern long dim1(void* v);
extern long dim2(void* v);
extern void femlib_ary1(void** v,long n_1, size_t o);
extern void femlib_ary2(void** v,long m_1, long n_1, size_t o);
extern double* lvec(void *A, long i, long j);

#endif
