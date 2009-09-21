#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ double x, y; char *label;} xyc;
typedef struct{ int a,b,c,A,B,C;} nde;

#define fp2mesh(fp,Zp,Np) estiva_fp2mesh(fp,Zp,Np)
#define plt(fp,Mid,Z,N,u)     estiva_plt(fp,Mid,Z,N,u)
#define np1(Z,N)          estiva_np1(Z,N)
#define forgamma(Z,i,str) for(estiva_initgamma();estiva_gamma(Z,&(i),str);)


extern void estiva_fp2mesh(FILE* fp, xyc** Zp, nde** Np);
extern void estiva_plt(FILE *fp, xyc *Mid, xyc *Z, nde *N, double *u);
extern void *estiva_np1(xyc *Z, nde *N);
extern void estiva_initgamma(void);
extern int  estiva_gamma(xyc* Z,int *p,char *str);


#endif
