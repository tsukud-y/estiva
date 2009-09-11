#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ double x, y; char *label;} xyc;
typedef struct{ int a,b,c,A,B,C;} nde;

#define fp2mesh(fp,Zp,Np) estiva_fp2mesh(fp,Zp,Np)
#define plt(fp,Z,N,u)     estiva_plt(fp,Z,N,u)

extern void estiva_fp2mesh(FILE* fp, xyc** Zp, nde** Np);
extern void estiva_plt(FILE *fp, xyc *Z, nde *N, double *u);

#endif
