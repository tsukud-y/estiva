#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ double x, y; char *label;} xyc;
typedef struct{ int a,b,c,A,B,C;} nde;

void estiva_fp2mesh(FILE* fp, xyc** Zp, nde** Np);
void estiva_plt(FILE *fp, xyc *Mid, xyc *Z, nde *N, double *u);
void *estiva_np1(xyc *Z, nde *N);
void estiva_initgamma(void);
int  estiva_gamma(xyc* Z,int *p,char *str);
xyc *estiva_fp2Z(FILE *fp);
void estiva_delaunay(xyc **Zo, nde **No);

#define fp2mesh(fp,Zp,Np) estiva_fp2mesh(fp,Zp,Np)
#define plt(fp,Mid,Z,N,u) estiva_plt(fp,Mid,Z,N,u)
#define np1(Z,N)          estiva_np1(Z,N)
#define forgamma(Z,i,str) for(estiva_initgamma();estiva_gamma(Z,&(i),str);)
#define fp2Z(fp)          estiva_fp2Z(fp)
#define delaunay(Z,N)     estiva_delaunay(&(Z),&(N)) 

#endif
