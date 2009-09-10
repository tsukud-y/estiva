#ifndef _MSH_H_
#define _MSH_H_

#include <stdio.h>

typedef struct{
  double x, y;
  char *label;
} xyc;

typedef struct{
  int a,b,c,A,B,C;
} nde;

extern void numgrid_initgamma(void);
extern int numgrid_gamma(xyc* Z,int *p,char *str);
#define forgamma(Z,i,str) for(numgrid_initgamma();numgrid_gamma(Z,&(i),str);)

extern void fp2msh(FILE* fp, xyc** Zp, nde** Np);
extern void *Ver2Mid(xyc *Z, nde *N);
extern double *S_(xyc *Z, nde *N);
#endif
