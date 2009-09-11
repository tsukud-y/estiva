#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ double x, y; char *label;} xyc;
typedef struct{ int a,b,c,A,B,C;} nde;


extern void fp2mesh(FILE* fp, xyc** Zp, nde** Np);


#endif
