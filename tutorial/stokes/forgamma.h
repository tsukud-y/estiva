#ifndef _MSH_H_
#define _MSH_H_

#define forgamma(Z,i,str) for(numgrid_initgamma();numgrid_gamma(Z,&(i),str);)

extern void numgrid_initgamma(void);
extern int numgrid_gamma(xyc* Z,int *p,char *str);


#endif
