#ifndef _DYNAMARY_H_
#define _DYNAMARY_H_
extern void *confary_ary1(void *A, int n, int size);
extern void *confary_ary2(void *A, int m, int n, int size);
extern int   confary_dim1(void *p);
extern int   confary_dim2(void *p);
extern int   confary_siz(void *p);
extern void  confary_reset(void *A);

#define ary1(b,n)   confary_ary1(&(b),  n,sizeof( *(b)))
#define ary2(A,m,n) \
confary_ary2(&(A),m,n,sizeof(**(A)));\
{int i;for(i=1;i<=dim2(A);i++) A[i]= &A[0][i*n];}
#define dim1(p)     confary_dim1(p) 
#define dim2(p)     confary_dim2(p) 
#define siz(p)      confary_siz(p)  
#define reset(A)    confary_reset(A) 
#endif
