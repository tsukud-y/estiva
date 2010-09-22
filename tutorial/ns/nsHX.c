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


static double hxij(long i, long j, double B1, double B2, double B3)
{
  double hxij, alphaBi, betaBi, gammaBi;

  alphaBi =     b[i]*B1 +     c[i]*B2;
  betaBi  = 2.0*e[i]*B1 +     d[i]*B2;
  gammaBi =     d[i]*B1 + 2.0*f[i]*B2;
  
  hxij =
    +  12.0*(alphaBi*ad[j])
    +   4.0*( betaBi*ad[j] + alphaBi*bd[j] + gammaBi*ad[j] + alphaBi*cd[j])
    +   1.0*( betaBi*cd[j] + gammaBi*bd[j])
    +   2.0*( betaBi*bd[j] + gammaBi*cd[j])
    ;
  return hxij/12.0;
}



void estiva_HX(MX **HXp, double *S, xyc *Z, nde *N)
{
  int I, J, i, j, e, a, b, c, A, B, C, n;
  double delta;
  static MX *HX;
  static double  B1, B2, B3;

  initmx(HX,dimp2(N)+1,21);

  n = dim1(N);
  for ( e = 1; e <= n; e++ ) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    delta = S[e];
    i = 0;
    foreach(I) &a, &b, &c, &A, &B, &C, end {
      ++i; j=0;
      foreach(J) &a, &b, &c, end {
	B1 = Z[b].y - Z[c].y; B1 /= 2.0*delta;
	B2 = Z[c].y - Z[a].y; B2 /= 2.0*delta;
	B3 = Z[a].y - Z[b].y; B3 /= 2.0*delta;
	mx(HX,I,J) = delta*hxij(i,++j, B1, B2, B3);
      }
    }
  }
  *HXp = HX;
}


