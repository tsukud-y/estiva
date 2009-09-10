#ifndef _NUMGRID_H_
#define _NUMGRID_H_
struct element { int a,b,c,ua,ub,uc,ma,mb,mc; double x,y,s;};
struct vertex { double x, y; char *l;};
struct midpoint { char *l;};
extern struct element  *numgrid_Ele;
extern struct vertex   *numgrid_Ver;
extern struct midpoint *numgrid_Mid;
extern void numgrid_initgamma(void);
extern int  numgrid_gamma(int *p, char *str); 
extern void numgrid_readfile(char *name, char *name1, char *name2);
extern void numgrid_lambda(float **A, int i);
extern int numgrid_e;
extern void FEM(void *v0,...);

#define Ele numgrid_Ele
#define Ver numgrid_Ver
#define Mid numgrid_Mid
#define forgamma(i,str) for(numgrid_initgamma();numgrid_gamma(&(i),str);)
#define lambda(A,i)                numgrid_lambda(A,i)
#define readfile(name,name1,name2) numgrid_readfile(name,name1,name2)
#define forFEM(e) for(e=(numgrid_e=1);e<=dim1(Ele);e= ++numgrid_e)
#endif
