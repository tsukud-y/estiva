#ifndef _ESTIVA_RF4_H_
#define _ESTIVA_RF4_H_

void  estiva_std_R4(void *x, void *y);
void *estiva_std_f4(void *x);
void  estiva_std_Rdestroy4(void *x);
void *estiva_std_f24(long size, void *x);

#define    R4(x,y)     estiva_std_R4(x,y)
#define Rdestroy4(x)   estiva_std_Rdestroy4(x)
#define static_new4(type,x)   Rnew4(x,type)
#define static_free4(x)       Rdestroy4(x)
#define static_bind4(type,x) (*(type*)estiva_std_f24(sizeof(type),x))

#endif 
