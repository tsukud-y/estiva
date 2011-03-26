#ifndef _ESTIVA_ARY_H_
#define _ESTIVA_ARY_H_

void    estiva_ary1(void** v,long n_1, size_t o);
void    estiva_ary2(void** v,long m_1, long n_1, size_t o);
size_t  estiva_dim0(void* v);
long    estiva_dim1(void* v);
long    estiva_dim2(void* v);

#define ary1(b,n_1)        estiva_ary1((void *)&b,n_1,sizeof(*b))
#define ary2(A,m_1,n_1)    estiva_ary2((void *)&A,m_1,n_1,sizeof(**A))
#define dim0(v)            estiva_dim0(v)
#define dim1(v)            estiva_dim1(v)
#define dim2(v)            estiva_dim2(v)

#endif
