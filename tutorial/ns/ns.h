#ifndef _ESTIVA_NS_H_
#define _ESTIVA_NS_H_

#include <stdio.h>
#include "estiva/mx.h"
#include "estiva/mesh.h"

#define nsM(M,S,N)            estiva_nsM(&M,S,N)
#define nsAX(AX,U,S,Z,N)      estiva_AX(&AX,U,S,Z,N)
#define nsAY(AY,V,S,Z,N)      estiva_AX(&AY,V,S,Z,N)
#define nsD(D,S,Z,N)          estiva_D(&D,S,Z,N)
#define nsHX(HX,S,Z,N)         estiva_HX(&HX,S,Z,N)
#define nsHY(HY,S,Z,N)         estiva_HY(&HY,S,Z,N)


void estiva_nsM(MX **Mp, double *S, nde *N);
void estiva_AX(MX **AXp, double *U, double *S, xyc *Z, nde *N);
void estiva_AY(MX **AYp, double *V, double *S, xyc *Z, nde *N);
void estiva_D(MX **Dp, double *S, xyc *Z, nde *N);
void estiva_HX(MX **HXp, double *S, xyc *Z, nde *N);
void estiva_HY(MX **HYp, double *S, xyc *Z, nde *N);
void setBCD(double b1, double b2, double c1, double c2, double s);
double alphaB(long j);
double  betaB(long j);
double gammaB(long j);
double alphaC(long j);
double  betaC(long j);
double gammaC(long j);
  
extern double Delta;

/*
//                                   a1    a2    a3    a4    a5    a6
static double a[] = { 0.0,          0.0,  0.0,  1,0,  0.0,  0.0,  0.0};
static double b[] = { 0.0,         -1.0,  0.0, -3.0,  0.0,  4.0,  0.0};
static double c[] = { 0.0,          0.0, -1.0, -3.0,  4.0,  0.0,  0.0};
static double d[] = { 0.0,          0.0,  0.0,  4.0, -4.0, -4.0,  4.0};
static double e[] = { 0.0,          2.0,  0.0,  2.0,  0.0, -4.0,  0.0};
static double f[] = { 0.0,          0.0,  2.0,  2.0, -4.0,  0.0,  0.0};


//                                  a'1   a'2   a'3
static double ad[] = { 0.0,         0.0,  0.0,  1.0};
static double bd[] = { 0.0,         1.0,  0.0, -1.0};
static double cd[] = { 0.0,         0.0,  1.0, -1.0};
*/

#endif
