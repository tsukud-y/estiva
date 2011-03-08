#ifndef _ESTIVA_TAYLORHOOD_H_
#define _ESTIVA_TAYLORHOOD_H_
#include <stdio.h>
#include "estiva/mx.h"
#include "estiva/mesh.h"

extern double estiva_TaylorHood_Delta;

double  estiva_TaylorHood_a(long i);
double  estiva_TaylorHood_b(long i);
double  estiva_TaylorHood_c(long i);
double  estiva_TaylorHood_d(long i);
double  estiva_TaylorHood_e(long i);
double  estiva_TaylorHood_f(long i);
double  estiva_TaylorHood_ad(long j);
double  estiva_TaylorHood_bd(long j);
double  estiva_TaylorHood_cd(long j);
double  estiva_TaylorHood_alphaB(long j);
double  estiva_TaylorHood_betaB(long j);
double  estiva_TaylorHood_gammaB(long j);
double  estiva_TaylorHood_alphaC(long j);
double  estiva_TaylorHood_betaC(long j);
double  estiva_TaylorHood_gammaC(long j);
void    estiva_TaylorHood_K(MX **K, long w);
void    estiva_TaylorHood_M(MX **M, long w);
void    estiva_TaylorHood_Hx(MX **Hx, long w);
void    estiva_TaylorHood_Hy(MX **Hy, long w);
void    estiva_genP2P2mx(MX **Mp, double (*func)(long i, long j), long w);
void    estiva_genP2P1mx(MX **Mp, double (*func)(long i, long j), long w);
void    estiva_setBCD(double b1, double b2, double c1, double c2, double s);  
double *estiva_setmesh(xyc *Z, nde *N);
void   estiva_TaylorHood_Ax(MX **Axp, double *U, long w);
void   estiva_TaylorHood_Ay(MX **Ayp, double *V, long w);
void   estiva_getZNS(xyc **Zp, nde **Np, double **S);

#define a(i)                  estiva_TaylorHood_a(i) 
#define b(i)                  estiva_TaylorHood_b(i) 
#define c(i)                  estiva_TaylorHood_c(i) 
#define d(i)                  estiva_TaylorHood_d(i) 
#define e(i)                  estiva_TaylorHood_e(i) 
#define f(i)                  estiva_TaylorHood_f(i) 
#define ad(j)                 estiva_TaylorHood_ad(j)
#define bd(j)                 estiva_TaylorHood_bd(j)
#define cd(j)                 estiva_TaylorHood_cd(j)
#define alphaB(j)             estiva_TaylorHood_alphaB(j)
#define  betaB(j)             estiva_TaylorHood_betaB(j)
#define gammaB(j)             estiva_TaylorHood_gammaB(j)
#define alphaC(j)             estiva_TaylorHood_alphaC(j)
#define  betaC(j)             estiva_TaylorHood_betaC(j)
#define gammaC(j)             estiva_TaylorHood_gammaC(j)
#define Delta                 estiva_TaylorHood_Delta
#define TaylorHood_K(K,w)     estiva_TaylorHood_K(&K,w)
#define TaylorHood_M(M,w)     estiva_TaylorHood_M(&M,w)
#define TaylorHood_Hx(M,w)    estiva_TaylorHood_Hx(&Hx,w)
#define TaylorHood_Hy(M,w)    estiva_TaylorHood_Hy(&Hy,w)
#define genP2P2mx(M,f,w)      estiva_genP2P2mx(&M,f,w)
#define genP2P1mx(M,f,w)      estiva_genP2P1mx(&M,f,w)
#define setBCD(b1,b2,c1,c2,s) estiva_setBCD(b1, b2, c1, c2, s) 
#define setmesh(Z,N)          estiva_setmesh(Z,N)
#define TaylorHood_Ax(Ax,U,w) estiva_TaylorHood_Ax(&Ax,U,w)
#define TaylorHood_Ay(Ay,V,w) estiva_TaylorHood_Ay(&Ay,V,w)
#define getZNS(Z,N,S)      estiva_getZNS(&Z,&N,&S)    

#endif
