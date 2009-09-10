#ifndef _ESTIVA_ARY_H_
#define _ESTIVA_ARY_H_

#include <stdlib.h>
#define dim0(v)         estiva_dim0(v)
#define dim1(v)         estiva_dim1(v)
#define dim2(v)         estiva_dim2(v)
#define ary1(b,n_1)     estiva_ary1((void **)&b,n_1,sizeof(*b))
#define ary2(A,m_1,n_1) estiva_ary2((void **)&A,m_1,n_1,sizeof(**A))

extern size_t estiva_dim0(void* v);
extern long   estiva_dim1(void* v);
extern long   estiva_dim2(void* v);
extern void   estiva_ary1(void** v,long n_1, size_t o);
extern void   estiva_ary2(void** v,long m_1, long n_1, size_t o);

#endif
