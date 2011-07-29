#ifndef _ESTIVA_RF5_H_
#define _ESTIVA_RF5_H_

void  estiva_std_R5(void *x, void *y);
void *estiva_std_f5(void *x);
void  estiva_std_Rdestroy5(void *x);
void *estiva_std_f25(long size, void *x);

#define    R5(x,y)     estiva_std_R5(x,y)
#define Rdestroy5(x)   estiva_std_Rdestroy5(x)
#define static_new5(type,x)   Rnew5(x,type)
#define static_free5(x)       Rdestroy5(x)
#define static_bind5(type,x) (*(type*)estiva_std_f25(sizeof(type),x))

#endif 
