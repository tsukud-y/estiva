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
#define a(i)                  estiva_a(i) 
#define b(i)                  estiva_b(i) 
#define c(i)                  estiva_c(i) 
#define d(i)                  estiva_d(i) 
#define e(i)                  estiva_e(i) 
#define f(i)                  estiva_f(i) 
#define ad(j)                 estiva_ad(j)
#define bd(j)                 estiva_bd(j)
#define cd(j)                 estiva_cd(j)
#define alphaB(j)             estiva_alphaB(j)
#define  betaB(j)             estiva_betaB(j)
#define gammaB(j)             estiva_gammaB(j)
#define alphaC(j)             estiva_alphaC(j)
#define  betaC(j)             estiva_betaC(j)
#define gammaC(j)             estiva_gammaC(j)


void   estiva_nsM(MX **Mp, double *S, nde *N);
void   estiva_AX(MX **AXp, double *U, double *S, xyc *Z, nde *N);
void   estiva_AY(MX **AYp, double *V, double *S, xyc *Z, nde *N);
void   estiva_D(MX **Dp, double *S, xyc *Z, nde *N);
void   estiva_HX(MX **HXp, double *S, xyc *Z, nde *N);
void   estiva_HY(MX **HYp, double *S, xyc *Z, nde *N);
double estiva_a(long i);
double estiva_b(long i);
double estiva_c(long i);
double estiva_d(long i);
double estiva_e(long i);
double estiva_f(long i);
double estiva_ad(long j);
double estiva_bd(long j);
double estiva_cd(long j);
double estiva_alphaB(long j);
double estiva_betaB(long j);
double estiva_gammaB(long j);
double estiva_alphaC(long j);
double estiva_betaC(long j);
double estiva_gammaC(long j);

void   setBCD(double b1, double b2, double c1, double c2, double s);  
extern double Delta;

#endif
