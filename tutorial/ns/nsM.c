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


static double mij(long i, long j)
{
  double mij =
    180.0*(a[i]*a[j])
    +  60.0*(a[i]*b[j] + b[i]*a[j] + a[i]*c[j] + c[i]*a[j])
    +  15.0*(a[i]*d[j] + d[i]*a[j] + b[i]*c[j] + c[i]*b[j])
    +  30.0*(b[i]*b[j] + c[i]*c[j] + a[i]*e[j] + e[i]*a[j] + a[i]*f[j] + f[i]*a[j])
    +  18.0*(b[i]*e[j] + e[i]*b[j] + c[i]*f[j] + f[i]*c[j])
    +   6.0*(b[i]*d[j] + d[i]*b[j] + c[i]*e[j] + e[i]*c[j] + b[i]*f[j] + f[i]*b[j] + c[i]*d[j] + d[i]*c[j])
    +   2.0*(d[i]*d[j] + e[i]*f[j] + f[i]*e[j])
    +   3.0*(d[i]*e[j] + e[i]*d[j] + d[i]*f[j] + f[i]*d[j])
    +  12.0*(e[i]*e[j] + f[i]*f[j]);

  return mij/180.0;
}


void estiva_nsM(MX **Mp, double *S, nde *N)
{
  int I, J, i, j, e, a, b, c, A, B, C, n;
  double delta;
  static MX *M;
  initmx(M,dimp2(N)+1,20);

  n = dim1(N);
  for ( e = 1; e <= n; e++ ) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    delta = S[e];
    i = 0;
    foreach(I) &a, &b, &c, &A, &B, &C, end {
      ++i; j=0;
      foreach(J) &a, &b, &c, &A, &B, &C, end {
	mx(M,I,J) = delta*mij(i,++j);
      }
    }
  }
  *Mp = M;
}
