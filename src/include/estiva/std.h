#ifndef _ESTIVA_STD_H_
#define _ESTIVA_STD_H_

void  estiva_swap(void *A, void *B, int size1, int size2);
void  estiva_cp(void *A, void *B, int size1, int size2);
void  estiva_std_R(void *x, void *y);
void *estiva_std_f(void *x);
void  estiva_std_Rdestroy(void *x);
void *estiva_std_f2(long size, void *x);

#define forall(m,i,n) for(i=m;i<=n;i++)
#define  max(x,y)     (x>y?(x):(y))
#define  min(x,y)     (x<y?(x):(y))
#define swap(a,b)     estiva_swap(&(a),&(b),sizeof(a),sizeof(b))
#define   cp(a,b)     estiva_cp(&(a),&(b),sizeof(a),sizeof(b))
#define    R(x,y)     estiva_std_R(x,y)
#define Rdestroy(x)   estiva_std_Rdestroy(x)
#define static_new(type,x)   Rnew(x,type)
#define static_free(x)       Rdestroy(x)
#define static_bind(type,x) (*(type*)estiva_std_f2(sizeof(type),x))

#endif 
