#ifndef _ESTIVA_H_
#define _ESTIVA_H_
extern void estiva_swap(void *A, void *B, int size1, int size2);
extern void estiva_cp(void *A, void *B, int size1, int size2);
extern void estiva_initecho(void);
extern int  estiva_echo(char *str);

#define forall(m,i,n) for(i=m;i<=n;i++)
#define more(x,y)     (x>y?(x):(y))
#define less(x,y)     (x<y?(x):(y))
#define absv(x)       ((x)>0.0?(x):-(x))
#define swap(a,b)     estiva_swap(&(a),&(b),sizeof(a),sizeof(b))
#define   cp(a,b)     estiva_cp(&(a),&(b),sizeof(a),sizeof(b))
#define echo(str)     for(estiva_initecho();estiva_echo(str);)
#endif 
