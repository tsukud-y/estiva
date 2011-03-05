#ifndef _ESTIVA_NS_H_
#define _ESTIVA_NS_H_

#include <stdio.h>
#include "estiva/mx.h"
#include "estiva/mesh.h"

#define nsM(M,S,N)            estiva_nsM(&M,S,N)
#define nsAX(AX,U,S,Z,N)      estiva_AX(&AX,U,S,Z,N)
#define nsAY(AY,V,S,Z,N)      estiva_AY(&AY,V,S,Z,N)
#define nsD(D,S,Z,N)          estiva_D(&D,S,Z,N)
#define nsHX(HX,S,Z,N)        estiva_HX(&HX,S,Z,N)
#define nsHY(HY,S,Z,N)        estiva_HY(&HY,S,Z,N)
#define boundary_condition(Z,N,A,b) estiva_boundary_condition(Z,N,A,b)
#define nsRhs(b,M,x)          estiva_nsRhs(b,M,x)
#define nsA(A,x,b,Z,N,K,M,Hx,Hy,AX,AY,t)  estiva_nsA(&A,x,b,Z,N,K,M,Hx,Hy,AX,AY,t)
#define genP2P2mx(M,f)        estiva_genP2P2mx(&M,f)
#define genP2P1mx(M,f)        estiva_genP2P1mx(&M,f)

void   estiva_nsM(MX **Mp, double *S, nde *N);
void   estiva_AX(MX **AXp, double *U, double *S, xyc *Z, nde *N);
void   estiva_AY(MX **AYp, double *V, double *S, xyc *Z, nde *N);
void   estiva_D(MX **Dp, double *S, xyc *Z, nde *N);
void   estiva_HX(MX **HXp, double *S, xyc *Z, nde *N);
void   estiva_HY(MX **HYp, double *S, xyc *Z, nde *N);
void   estiva_boundary_condition(xyc *Z, nde *N, MX *A, double *b);
void   estiva_nsRhs(double *b, MX *M, double *x);
void   estiva_nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N, MX *K, MX *M, MX *Hx, MX *Hy, MX *AX, MX *AY, double tau);
void   setBCD(double b1, double b2, double c1, double c2, double s);  

void estiva_genP2P2mx(MX **Mp, double (*func)(long i, long j));
void estiva_genP2P1mx(MX **Mp, double (*func)(long i, long j));

double  mij(long i, long j);
double  kij(long i, long j);
double hxij(long i, long j);
double hyij(long i, long j);

#endif
