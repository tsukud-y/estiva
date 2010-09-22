#include "ns.h"
#include "fem.h"


#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/mx.h"

//                                   a1    a2    a3    a4    a5    a6
static double b[] = { 0.0,         -1.0,  0.0, -3.0,  0.0,  4.0,  0.0};
static double c[] = { 0.0,          0.0, -1.0, -3.0,  4.0,  0.0,  0.0};
static double d[] = { 0.0,          0.0,  0.0,  4.0, -4.0, -4.0,  4.0};
static double e[] = { 0.0,          2.0,  0.0,  2.0,  0.0, -4.0,  0.0};
static double f[] = { 0.0,          0.0,  2.0,  2.0, -4.0,  0.0,  0.0};

//                                  a'1   a'2   a'3
static double ad[] = { 0.0,         0.0,  0.0,  1.0};
static double bd[] = { 0.0,         1.0,  0.0, -1.0};
static double cd[] = { 0.0,         0.0,  1.0, -1.0};


static double hyij(long i, long j, double C1, double C2, double C3)
{
  double hyij, alphaCi, betaCi, gammaCi;

  alphaCi =     b[i]*C1 +     c[i]*C2;
  betaCi  = 2.0*e[i]*C1 +     d[i]*C2;
  gammaCi =     d[i]*C1 + 2.0*f[i]*C2;
  
  hyij =
    +  12.0*(alphaCi*ad[j])
    +   4.0*( betaCi*ad[j] + alphaCi*bd[j] + gammaCi*ad[j] + alphaCi*cd[j])
    +   1.0*( betaCi*cd[j] + gammaCi*bd[j])
    +   2.0*( betaCi*bd[j] + gammaCi*cd[j])
    ;
  return hyij/12.0;
}


void estiva_HY(MX **HYp, double *S, xyc *Z, nde *N)
{
  int I, J, i, j, e, a, b, c, A, B, C, n;
  double delta;
  static MX *HY;
  static double  C1, C2, C3;

  initmx(HY,dimp2(N)+1,21);

  n = dim1(N);
  for ( e = 1; e <= n; e++ ) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    delta = S[e];
    i = 0;
    foreach(I) &a, &b, &c, &A, &B, &C, end {
      ++i; j=0;
      foreach(J) &a, &b, &c, end {
	C1 = Z[b].x - Z[c].x; C1 /= 2.0*delta;
	C2 = Z[c].x - Z[a].x; C2 /= 2.0*delta;
	C3 = Z[a].x - Z[b].x; C3 /= 2.0*delta;
	mx(HY,I,J) = delta*hyij(i,++j, C1, C2, C3);
      }
    }
  }
  *HYp = HY;
}


