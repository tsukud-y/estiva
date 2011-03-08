#ifndef _ESTIVA_NS_H_
#define _ESTIVA_NS_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/mx.h"
#include "estiva/mesh.h"
#include "estiva/TaylorHood.h"

void   estiva_boundary_condition(MX *A, double *b);
void   estiva_nsRhs(double *b, MX *M, double *x);
void   estiva_nsA(MX **Ap, double *x, double *b, MX *K, MX *M, MX *Hx, MX *Hy, MX *AX, MX *AY, double tau, long w );

#define boundary_condition(A,b)          estiva_boundary_condition(A,b)
#define nsRhs(b,M,x)                     estiva_nsRhs(b,M,x)
#define nsA(A,x,b,K,M,Hx,Hy,Ax,Ay,tau,w) estiva_nsA(&A,x,b,K,M,Hx,Hy,Ax,Ay,tau,w)
#endif
