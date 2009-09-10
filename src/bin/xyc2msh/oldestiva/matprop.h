#ifndef _MATPROP_H_
#define _MATPROP_H_
extern int matprop_halfbw(void *A);

#define halfbw(A) matprop_halfbw(A)
#endif 
