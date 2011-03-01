#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ double x, y; char *label;} xyc;
typedef struct{ int a,b,c,A,B,C;} nde;

void estiva_fp2xyc(FILE *fp, xyc **Zp);
void estiva_fp2mesh(FILE* fp, xyc **Zp, nde **Np);
void estiva_plt(FILE *fp, xyc *Mid, xyc *Z, nde *N, double *u);
void *estiva_np1(xyc *Z, nde *N);
void estiva_initgamma(void);
int  estiva_gamma(xyc* Z,int *p,char *str);
void estiva_delaunay(xyc **Zo, nde **No);
void estiva_fprintmesh(FILE *fp, xyc *Z, nde *N);
void estiva_xmesh(xyc *Z);
void estiva_pushxyc(void *qp, double x, double y, char *label);
void estiva_genmesh(void *qp, xyc **Zp, nde **Np);
void estiva_p2(xyc *Z, nde *N);
void estiva_rectmesh(xyc **Zp, nde **Np);
long estiva_dimp2(nde *N);
void estiva_pltp2(double *x, xyc * Z, nde *N);
void estiva_forgammap2(xyc *Z, nde *N, char *label);
int estiva_forgammap2_loop(long *ip);

#define fp2xyc(fp,Z)          estiva_fp2xyc(fp,&(Z))
#define fp2mesh(fp,Z,N)       estiva_fp2mesh(fp,&(Z),&(N))
#define plt(fp,Mid,Z,N,u)     estiva_plt(fp,Mid,Z,N,u)
#define np1(Z,N)              estiva_np1(Z,N)
#define forgamma(Z,i,str)     for(estiva_initgamma();estiva_gamma(Z,&(i),str);)
#define delaunay(Z,N)         estiva_delaunay(&(Z),&(N)) 
#define fprintmesh(fp,Z,N)    estiva_fprintmesh(fp,Z,N)
#define xmesh(Z)              estiva_xmesh(Z)
#define pushxyc(q,x,y,label)  estiva_pushxyc(&(q),x,y,label)
#define genmesh(q,Z,N)        estiva_genmesh(&(q),&(Z),&(N))
#define p2(Z,N)               estiva_p2(Z,N)
#define rectmesh(Z,N)         estiva_rectmesh(&Z,&N)
#define dimp2(N)              estiva_dimp2(N)
#define pltp2(x,Z,N)          estiva_pltp2(x, Z, N)
#define forgammap2(i,label,Z,N)						\
  for ( estiva_forgammap2(Z,N,label); estiva_forgammap2_loop(&(i));)

#endif
