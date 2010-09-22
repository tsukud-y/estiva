#include "ns.h"
#include "fem.h"


#include <stdio.h>
#include <estiva/mx.h>
#include <estiva/mesh.h>
#include <estiva/foreach.h>
#include <estiva/ary.h>

//                                   a1    a2    a3    a4    a5    a6
static double b[] = { 0.0,         -1.0,  0.0, -3.0,  0.0,  4.0,  0.0};
static double c[] = { 0.0,          0.0, -1.0, -3.0,  4.0,  0.0,  0.0};
static double d[] = { 0.0,          0.0,  0.0,  4.0, -4.0, -4.0,  4.0};
static double e[] = { 0.0,          2.0,  0.0,  2.0,  0.0, -4.0,  0.0};
static double f[] = { 0.0,          0.0,  2.0,  2.0, -4.0,  0.0,  0.0};


static double dij(long i, long j, double B1, double B2, double B3, double C1, double C2, double C3)
{
  double dij, alphaBi, betaBi, gammaBi, alphaCi, betaCi, gammaCi;
  double      alphaBj, betaBj, gammaBj, alphaCj, betaCj, gammaCj;

  alphaBi =     b[i]*B1 +     c[i]*B2;
  betaBi  = 2.0*e[i]*B1 +     d[i]*B2;
  gammaBi =     d[i]*B1 + 2.0*f[i]*B2;

  alphaBj =     b[j]*B1 +     c[j]*B2;
  betaBj  = 2.0*e[j]*B1 +     d[j]*B2;
  gammaBj =     d[j]*B1 + 2.0*f[j]*B2;
  
  alphaCi =     b[i]*C1 +     c[i]*C2;
  betaCi  = 2.0*e[i]*C1 +     d[i]*C2;
  gammaCi =     d[i]*C1 + 2.0*f[i]*C2;
  
  alphaCj =     b[j]*C1 +     c[j]*C2;
  betaCj  = 2.0*e[j]*C1 +     d[j]*C2;
  gammaCj =     d[j]*C1 + 2.0*f[j]*C2;

  dij =
    +  12.0*(alphaBi*alphaBj + alphaCi*alphaCj)
    +   4.0*(alphaBi* betaBj +  betaBi*alphaBj + alphaBi*gammaBj * gammaBi*alphaBj)
    +   4.0*(alphaCi* betaCj +  betaCi*alphaCj + alphaCi*gammaCj * gammaCi*alphaCj)
    +   1.0*( betaBi*gammaBj + gammaBi* betaBj +  betaCi*gammaCj * gammaCi* betaCj)
    +   2.0*( betaBi* betaBj + gammaBi*gammaBj +  betaCi* betaCj * gammaCi*gammaCj)
    ;
  return dij/12.0;
}


void estiva_D(MX **Dp, double *S, xyc *Z, nde *N)
{
  int I, J, i, j, e, a, b, c, A, B, C, n;
  double delta;
  static MX *D;
  static double B1, B2, B3, C1, C2, C3;

  initmx(D,dimp2(N)+1,28);

  n = dim1(N);
  for ( e = 1; e <= n; e++ ) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    delta = S[e];
    i = 0;
    foreach(I) &a, &b, &c, &A, &B, &C, end {
      ++i; j=0;
      foreach(J) &a, &b, &c, &A, &B, &C, end {
        B1 = Z[b].y - Z[c].y; B1 /= 2.0*delta;
        B2 = Z[c].y - Z[a].y; B2 /= 2.0*delta;
        B3 = Z[a].y - Z[b].y; B3 /= 2.0*delta;
        C1 = Z[b].x - Z[c].x; C1 /= 2.0*delta;
        C2 = Z[c].x - Z[a].x; C2 /= 2.0*delta;
        C3 = Z[a].x - Z[b].x; C3 /= 2.0*delta;
        mx(D,I,J) = delta*dij(i,++j, B1, B2, B3, C1, C2, C3);
      }
    }
  }
  *Dp = D;
}
