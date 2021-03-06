#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/TaylorHood.h"


static double ayij(long i, long j, double *v, double C1, double C2, double C3)
{
  double ayij, Av, Bv, Cv, Dv, Ev, Fv, alphaCj, betaCj, gammaCj;
  long k;

  for (Av = 0.0, k = 1; k <= 6; k++) Av += a(k)*v[k];
  for (Bv = 0.0, k = 1; k <= 6; k++) Bv += b(k)*v[k];
  for (Cv = 0.0, k = 1; k <= 6; k++) Cv += c(k)*v[k];
  for (Dv = 0.0, k = 1; k <= 6; k++) Dv += d(k)*v[k];
  for (Ev = 0.0, k = 1; k <= 6; k++) Ev += e(k)*v[k];
  for (Fv = 0.0, k = 1; k <= 6; k++) Fv += f(k)*v[k];

  alphaCj =     b(j)*C1 +     c(j)*C2;
  betaCj  = 2.0*e(j)*C1 +     d(j)*C2;
  gammaCj =     d(j)*C1 + 2.0*f(j)*C2;

  ayij =
    + 420.0 * (a(i)*Av                              ) * (3.0*alphaCj +     betaCj +     gammaCj)
    + 105.0 * (a(i)*Bv + b(i)*Av                    ) * (4.0*alphaCj + 2.0*betaCj +     gammaCj)
    + 105.0 * (a(i)*Cv + c(i)*Av                    ) * (4.0*alphaCj +     betaCj + 2.0*gammaCj)
    +  21.0 * (a(i)*Dv + d(i)*Av + b(i)*Cv + c(i)*Bv) * (5.0*alphaCj + 2.0*betaCj + 2.0*gammaCj)
    +  42.0 * (a(i)*Ev + e(i)*Av + b(i)*Bv          ) * (5.0*alphaCj + 3.0*betaCj +     gammaCj)
    +  42.0 * (a(i)*Fv + f(i)*Av + c(i)*Cv          ) * (5.0*alphaCj +     betaCj + 3.0*gammaCj)
    +   7.0 * (b(i)*Dv + d(i)*Bv + c(i)*Ev + e(i)*Cv) * (6.0*alphaCj + 3.0*betaCj + 2.0*gammaCj)
    +   7.0 * (b(i)*Fv + f(i)*Bv + c(i)*Dv + d(i)*Cv) * (6.0*alphaCj + 2.0*betaCj + 3.0*gammaCj)
    +  21.0 * (b(i)*Ev + e(i)*Bv                    ) * (6.0*alphaCj + 4.0*betaCj +     gammaCj)
    +  21.0 * (c(i)*Fv + f(i)*Cv                    ) * (6.0*alphaCj +     betaCj + 4.0*gammaCj)
    +   3.0 * (d(i)*Ev + e(i)*Dv                    ) * (7.0*alphaCj + 4.0*betaCj + 2.0*gammaCj)
    +   3.0 * (d(i)*Fv + f(i)*Dv                    ) * (7.0*alphaCj + 2.0*betaCj + 4.0*gammaCj)
    +   2.0 * (e(i)*Fv + f(i)*Ev + d(i)*Dv          ) * (7.0*alphaCj + 3.0*betaCj + 3.0*gammaCj)
    +  12.0 * (e(i)*Ev                              ) * (7.0*alphaCj + 5.0*betaCj +     gammaCj)
    +  12.0 * (f(i)*Fv                              ) * (7.0*alphaCj +     betaCj + 5.0*gammaCj)
    ;
  return ayij/1260.0;
}


void estiva_TaylorHood_Ay(MX **AYp, double *V, long w)
{
  int I, J, i, j, e, a, b, c, A, B, C, n;
  double delta;
  static MX *AY;
  static double *v, C1, C2, C3;
  static xyc *Z; static nde *N; static double *S;

  getZNS(Z,N,S);
  ary1(v,7);

  initmx(AY,dimp2(N)+1,w);
  clearmx(AY);

  n = dim1(N);
  for ( e = 1; e <= n; e++ ) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    delta = S[e];
    i = 0;
    foreach(I, &a, &b, &c, &A, &B, &C) {
      ++i; j=0;
      foreach(J, &a, &b, &c, &A, &B, &C) {
        C1 = Z[c].x - Z[b].x; C1 /= 2.0*delta;
        C2 = Z[a].x - Z[c].x; C2 /= 2.0*delta;
        C3 = Z[b].x - Z[a].x; C3 /= 2.0*delta;
        v[1] = V[a];
        v[2] = V[b];
        v[3] = V[c];
        v[4] = V[A];
        v[5] = V[B];
        v[6] = V[C];
        mx(AY,I,J) += delta*ayij(i,++j, v, C1, C2, C3);
      }
    }
  }
  *AYp = AY;
}
