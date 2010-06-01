#ifndef _ESTIVA_STD_H_
#define _ESTIVA_STD_H_
void estiva_swap(void *A, void *B, int size1, int size2);
void estiva_cp(void *A, void *B, int size1, int size2);

#define forall(m,i,n) for(i=m;i<=n;i++)
#define  max(x,y)     (x>y?(x):(y))
#define  min(x,y)     (x<y?(x):(y))
#define swap(a,b)     estiva_swap(&(a),&(b),sizeof(a),sizeof(b))
#define   cp(a,b)     estiva_cp(&(a),&(b),sizeof(a),sizeof(b))
#endif 
