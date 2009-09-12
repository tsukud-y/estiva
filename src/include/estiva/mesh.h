#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ double x, y; char *label;} xyc;
typedef struct{ int a,b,c,A,B,C;} nde;

#define fp2mesh(fp,Zp,Np) estiva_fp2mesh(fp,Zp,Np)
#define plt(fp,Z,N,u)     estiva_plt(fp,Z,N,u)
#define pltuv(fp,Z,u,v)     estiva_pltuv(fp,Z,u,v)
#define forgamma(Z,i,str) for(estiva_initgamma();estiva_gamma(Z,&(i),str);)

extern void estiva_fp2mesh(FILE* fp, xyc** Zp, nde** Np);
extern void estiva_plt(FILE *fp, xyc *Z, nde *N, double *u);
extern void estiva_pltuv(FILE *fp, xyc *Z, double *u, double *v);
extern void estiva_initgamma(void);
extern int  estiva_gamma(xyc* Z,int *p,char *str);


#endif
