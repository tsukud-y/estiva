#ifndef RF_H
#define RF_H
extern void estiva_R(void *xi, void *yi);
#define  R(x,y) estiva_R(x,y) 
extern void *estiva_f(void *xi);
#define  f(x)   estiva_f(x) 
#endif
