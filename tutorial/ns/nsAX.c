#include "ns.h"
#include "fem.h"


#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/mx.h"

//                                   a1    a2    a3    a4    a5    a6
static double a[] = { 0.0,          0.0,  0.0,  1,0,  0.0,  0.0,  0.0};
static double b[] = { 0.0,         -1.0,  0.0, -3.0,  0.0,  4.0,  0.0};
static double c[] = { 0.0,          0.0, -1.0, -3.0,  4.0,  0.0,  0.0};
static double d[] = { 0.0,          0.0,  0.0,  4.0, -4.0, -4.0,  4.0};
static double e[] = { 0.0,          2.0,  0.0,  2.0,  0.0, -4.0,  0.0};
static double f[] = { 0.0,          0.0,  2.0,  2.0, -4.0,  0.0,  0.0};


static double axij(long i, long j, double *u, double B1, double B2, double B3)
{
  double axij, Au, Bu, Cu, Du, Eu, Fu, alphaBj, betaBj, gammaBj;
  long k;

  for (Au = 0.0, k = 1; k <= 6; k++) Au += a[k]*u[k];
  for (Bu = 0.0, k = 1; k <= 6; k++) Bu += b[k]*u[k];
  for (Cu = 0.0, k = 1; k <= 6; k++) Cu += c[k]*u[k];
  for (Du = 0.0, k = 1; k <= 6; k++) Du += d[k]*u[k];
  for (Eu = 0.0, k = 1; k <= 6; k++) Eu += e[k]*u[k];
  for (Fu = 0.0, k = 1; k <= 6; k++) Fu += f[k]*u[k];

  alphaBj =     b[j]*B1 +     c[j]*B2;
  betaBj  = 2.0*e[j]*B1 +     d[j]*B2;
  gammaBj =     d[j]*B1 + 2.0*f[j]*B2;
  
  axij =
    + 420.0 * (a[i]*Au                              ) * (3.0*alphaBj +     betaBj +     gammaBj)
    + 105.0 * (a[i]*Bu + b[i]*Au                    ) * (4.0*alphaBj + 2.0*betaBj +     gammaBj)
    + 105.0 * (a[i]*Cu + c[i]*Au                    ) * (4.0*alphaBj +     betaBj + 2.0*gammaBj)
    +  21.0 * (a[i]*Du + d[i]*Au + b[i]*Cu + c[i]*Bu) * (5.0*alphaBj + 2.0*betaBj + 2.0*gammaBj)
    +  42.0 * (a[i]*Eu + e[i]*Au + b[i]*Bu          ) * (5.0*alphaBj + 3.0*betaBj +     gammaBj)
    +  42.0 * (a[i]*Fu + f[i]*Au + c[i]*Cu          ) * (5.0*alphaBj +     betaBj + 3.0*gammaBj)
    +   7.0 * (b[i]*Du + d[i]*Bu + c[i]*Eu + e[i]*Cu) * (6.0*alphaBj + 3.0*betaBj + 2.0*gammaBj)
    +   7.0 * (b[i]*Fu + f[i]*Bu + c[i]*Du + d[i]*Cu) * (6.0*alphaBj + 2.0*betaBj + 3.0*gammaBj)
    +  21.0 * (b[i]*Eu + e[i]*Bu                    ) * (6.0*alphaBj + 4.0*betaBj +     gammaBj)
    +  21.0 * (c[i]*Fu + f[i]*Cu                    ) * (6.0*alphaBj +     betaBj + 4.0*gammaBj)
    +   3.0 * (d[i]*Eu + e[i]*Du                    ) * (7.0*alphaBj + 4.0*betaBj + 2.0*gammaBj)
    +   3.0 * (d[i]*Fu + f[i]*Du                    ) * (7.0*alphaBj + 2.0*betaBj + 4.0*gammaBj)
    +   2.0 * (e[i]*Fu + f[i]*Eu + d[i]*Du          ) * (7.0*alphaBj + 3.0*betaBj + 3.0*gammaBj)
    +  12.0 * (e[i]*Eu                              ) * (7.0*alphaBj + 5.0*betaBj +     gammaBj)
    +  12.0 * (f[i]*Fu                              ) * (7.0*alphaBj +     betaBj + 5.0*gammaBj)
    ;
  return axij/1260.0;
}

void estiva_AX(MX **AXp, double *U, double *S, xyc *Z, nde *N)
{
  int I, J, i, j, e, a, b, c, A, B, C, n;
  double delta;
  static MX *AX;
  static double *u, B1, B2, B3;

  ary1(u,7);

  initmx(AX,dimp2(N)+1,28);

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
	u[1] = U[a];
	u[2] = U[b];
	u[3] = U[c];
	u[4] = U[A];
	u[5] = U[B];
	u[6] = U[C];
	mx(AX,I,J) = delta*axij(i,++j, u, B1, B2, B3);
      }
    }
  }
  *AXp = AX;
}


