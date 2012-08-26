#ifndef _estivaplus_h
#define _estivaplus_h

#define NDEBUG
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>

typedef struct{ int a,b,c,A,B,C;} nde;
typedef struct{ double x, y; char *label;} xyc;
typedef struct{double **A;long **IA;long n,w,I,J;double a;} MX;
typedef struct { void *elem; void *next;} que;

typedef std::vector< std::map<unsigned int, double> > Matrix;
typedef std::vector<double> Vector;
extern  std::map<unsigned int, double>::iterator estiva_aiterator;

void   PlotP2(Vector &xv);
void   Matrix2mx(Matrix &a, MX **Ap);
void   Vector2ary(std::vector<double> &bv, double *b);
void   nsRhs(Vector &b, Matrix &M, Vector &x);
long   getMatrixsize(Matrix &Hxm);
Matrix mx2Matrix(MX *A);
Matrix nsA(Matrix &M,Matrix &K,Matrix &Ax,Matrix &Ay,Matrix &Hx,Matrix &Hey);
Vector ary2Vector(double *x);
double tau();
Matrix TaylorHood_M();
Matrix TaylorHood_K();
Matrix TaylorHood_Hx();
Matrix TaylorHood_Hy();
Matrix TaylorHood_Ax(Vector &xv);
Matrix TaylorHood_Ay(Vector &xv);
void   BoundaryCondition(Matrix &Am, Vector &bv);
void   Solver(Matrix &A, Vector &x, Vector &b);
void   MatrixClear(Matrix &A);
void   MatrixDisp(Matrix &A);
void   estiva_putnode(double x, double y, char *label);
void   GenerateMesh(void);
void   Putnode(double x, double y, std::string label);
void   GenMeshP(void);
void PoiseuilleMesh(double h, double W, double H);

extern "C" {
  int     estiva_initop(int*,char***);
  void    cylindermesh(xyc**,nde**);
  void    estiva_fprintmesh(FILE*,xyc*,nde*);
  long    estiva_dimp2(void*);
  long    estiva_dim1(void*);
  void    estiva_ary1(void**,long,size_t);
  void    estiva_setmesh(xyc*,nde*);
  void    estiva_TaylorHood_M(MX **,long);
  void    estiva_TaylorHood_K(MX **,long);
  void    estiva_TaylorHood_Hx(MX **,long);
  void    estiva_TaylorHood_Hy(MX **,long);
  void    estiva_TaylorHood_Ax(MX **,double *,long);
  void    estiva_TaylorHood_Ay(MX **,double *,long);
  void    estiva_pltmx(MX *);
  int     estiva_defop(char*);
  char   *estiva_getop(char*);
  void    estiva_nsRhs(double*,MX*,double*);
  void    estiva_qmrsolver(MX*,double*,double*);
  void    estiva_cgssolver(MX*,double*,double*);
  void    estiva_bicgstabsolver(MX*,double*,double*);
  void    estiva_pltp2(double*,xyc*,nde*);
  void    estiva_fprintvec(double*,char*);
  void    estiva_boundary_condition(MX*,double*);
  double *estiva_mx(MX*,long,long);
  void    estiva_getZNS(xyc**,nde**,double**);
  void    estiva_initmx(MX**,long,long);
  void    estiva_clearmx(MX*);
  void    estiva_fornonzeromx(MX*);
  int     estiva_fornonzeromx_loop(MX *A, long *Ip, long *Jp);
  double *estiva_mx(MX*,long,long);
  void    estiva_forgammap1(long *ip);
  int     estiva_forgammap1_loop(long *ip, char *NAME, xyc *Z);
  void    estiva_forgammap2(long *ip, xyc *Z, nde *N, char *label);
  int     estiva_forgammap2_loop(long *ip);
  void    estiva_zerofillrow(MX *A, long i);
  void    estiva_fprintmesh(FILE *, xyc *, nde *);
  void  estiva_enq(  que  *q, void   *e, size_t n);
  void  estiva_initq(que **q);
  void *estiva_lup(  que  *q, void **s2, size_t n);
  void  estiva_push( que **q, void   *e, size_t n);
  void *estiva_pop(  que **q, void   *e, size_t n);
  void  estiva_forq(que *q, void **e);
  int   estiva_forq_loop(void **e);
  void  estiva_pushxyc(void *qp, double x, double y, char *label);
  void  estiva_genmesh(void *qp, xyc **Zp, nde **Np);
  void  estiva_p2(xyc *Z, nde *N);
}

#define forMatrixNonzero(a)                                             \
  for (unsigned int i=0;estiva_aiterator=a[i].begin(),i<a.capacity();i++) \
    for (unsigned int j=0;						\
	 j = (*estiva_aiterator).first, estiva_aiterator != a[i].end();	\
	 estiva_aiterator++)

#define fornonzeromx(A,i,j)                                             \
  for ( estiva_fornonzeromx(A); estiva_fornonzeromx_loop(A,&(i),&(j));)

#define forgammap1(i,NAME,Z)                                            \
  for ( estiva_forgammap1(&(i)); estiva_forgammap1_loop(&(i),(char*)NAME,Z);)

#define forgammap2(i,label,Z,N)                                         \
  for ( estiva_forgammap2(&(i),Z,N,(char*)label);estiva_forgammap2_loop(&(i));)


#define dim1(x)           estiva_dim1(x)
#define setmesh(Z,N)      estiva_setmesh(Z,N)
#define initop(argc,argv) estiva_initop(&(argc),&(argv))
#define dimp2(N)          estiva_dimp2(N)
#define ary1(x,n)         estiva_ary1((void**)&x,n,sizeof(*x))
#define defop(f)          estiva_defop((char*)f)
#define getop(f)          estiva_getop((char*)f)
#define opi(kn,n)         (kn = defop("-"#kn) ? atoi(getop("-"#kn)) : n)
#define opf(kn,n)         (kn = defop("-"#kn) ? atof(getop("-"#kn)) : n)
#define mx(A,i,j)                  (*estiva_mx(A,i,j))
#define zerofillrow(A,i)             estiva_zerofillrow(A,i)
#define fprintmesh(fp,Z,N)  estiva_fprintmesh(fp,Z,N)
#define enq(q,e)      estiva_enq(   q, (void*)&e, sizeof(e))
#define initq(q)      estiva_initq(&q)
#define lup(q,e)      estiva_lup(   q, (void*)&e, sizeof(e))
#define pop(q,e)      estiva_pop(  &q, (void*)&e, sizeof(e))
#define push(q,e)     estiva_push( &q, (void*)&e, sizeof(e))
#define pushxyc(q,x,y,label)          estiva_pushxyc(&(q),x,y,(char*)label)
#define genmesh(q,Z,N)                estiva_genmesh(&(q),&(Z),&(N))
#define p2(Z,N)                       estiva_p2(Z,N)
#define putnode(x,y,label) estiva_putnode(x,y,(char*)label)
#define forq(q,e)                                                       \
  for( estiva_forq(q,(void*)&e); estiva_forq_loop((void*)&e);)

#endif
