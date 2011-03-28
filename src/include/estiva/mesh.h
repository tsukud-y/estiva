#ifndef _ESTIVA_MESH_H_
#define _ESTIVA_MESH_H_

typedef struct{ int a,b,c,A,B,C;} nde;
typedef struct{ double x, y; char *label;} xyc;

void    estiva_delaunay(xyc **Zo, nde **No);
long    estiva_dimp2(nde *N);
void    estiva_femdelta(double **Sp, xyc *Z, nde *N);
void    estiva_fp2xyc(void *fp, xyc **Zp);
void    estiva_fp2mesh(void* fp, xyc **Zp, nde **Np);
void    estiva_fprintmesh(void *fp, xyc *Z, nde *N);
void    estiva_genmesh(void *qp, xyc **Zp, nde **Np);
void   *estiva_np1(xyc *Z, nde *N);
void    estiva_p2(xyc *Z, nde *N);
void    estiva_plt(void *fp, xyc *Mid, xyc *Z, nde *N, double *u);
void    estiva_pltp2(double *x, xyc * Z, nde *N);
void    estiva_pushxyc(void *qp, double x, double y, char *label);
void    estiva_rectmesh(xyc **Zp, nde **Np);
void    estiva_xmesh(xyc *Z);
void    estiva_forgamma(void);
int     estiva_forgamma_loop(xyc* Z,int *p,char *str);
void    estiva_forgammap1(long *ip);
int     estiva_forgammap1_loop(long *ip, char *NAME, xyc *Z);
void    estiva_forgammap2(long *ip, xyc *Z, nde *N, char *label);
int     estiva_forgammap2_loop(long *ip);

#define delaunay(Z,N)                 estiva_delaunay(&(Z),&(N)) 
#define dimp2(N)                      estiva_dimp2(N)
#define femdelta(S,Z,N)               estiva_femdelta(&S,Z,N)
#define fp2xyc(fp,Z)                  estiva_fp2xyc(fp,&(Z))
#define fp2mesh(fp,Z,N)               estiva_fp2mesh(fp,&(Z),&(N))
#define fprintmesh(fp,Z,N)            estiva_fprintmesh(fp,Z,N)
#define genmesh(q,Z,N)                estiva_genmesh(&(q),&(Z),&(N))
#define np1(Z,N)                      estiva_np1(Z,N)
#define p2(Z,N)                       estiva_p2(Z,N)
#define plt(fp,Mid,Z,N,u)             estiva_plt(fp,Mid,Z,N,u)
#define pltp2(x,Z,N)                  estiva_pltp2(x, Z, N)
#define pushxyc(q,x,y,label)          estiva_pushxyc(&(q),x,y,label)
#define rectmesh(Z,N)                 estiva_rectmesh(&Z,&N)
#define xmesh(Z)                      estiva_xmesh(Z)
#define forgamma(Z,i,str)						\
  for ( estiva_forgamma();            estiva_forgamma_loop(Z,&(i),str);)
#define forgammap1(i,NAME,Z)						\
  for ( estiva_forgammap1(&(i));      estiva_forgammap1_loop(&(i),NAME,Z);)
#define forgammap2(i,label,Z,N)						\
  for ( estiva_forgammap2(&(i),Z,N,label); estiva_forgammap2_loop(&(i));)

#endif
