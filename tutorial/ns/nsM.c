#include "ns.h"
#include "fem.h"


#include <math.h>
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/mx.h"

static double a[] = { NAN,          0.0,  0.0,  1,0,  0.0,  0.0,  0.0, NAN};
static double b[] = { NAN,         -1.0,  0.0, -3.0,  0.0,  4.0,  0.0, NAN};
static double c[] = { NAN,          0.0, -1.0, -3.0,  4.0,  0.0,  0.0, NAN};
static double d[] = { NAN,          0.0,  0.0,  4.0, -4.0, -4.0,  4.0, NAN};
static double e[] = { NAN,          2.0,  0.0,  2.0,  0.0, -4.0,  0.0, NAN};
static double f[] = { NAN,          0.0,  2.0,  2.0, -4.0,  0.0,  0.0, NAN};


static double l(long i, long j)
{
  return (Delta/180.0)*
    (180.0*(a[i]*a[j])                                                                                     +
     60.0 *(a[i]*b[j] + b[i]*a[j] + a[i]*c[j] + c[i]*a[j])                                                 +
     15.0 *(a[i]*d[j] + d[i]*a[j] + b[i]*c[j] + c[i]*b[j])                                                 +
     30.0 *(b[i]*b[j] + c[i]*c[j] + a[i]*e[j] + e[i]*a[j] + a[i]*f[j] + f[i]*a[j])                         +
     18.0 *(b[i]*e[j] + e[i]*b[j] + c[i]*f[j] + f[i]*c[j])                                                 +
     6.0  *(b[i]*d[j] + d[i]*b[j] + c[i]*e[j] + e[i]*c[j] + b[i]*f[j] + f[i]*b[j] + c[i]*d[j] + d[i]*c[j]) +
     2.0  *(d[i]*d[j] + e[i]*f[j] + f[i]*e[j])                                                             +
     3.0  *(d[i]*e[j] + e[i]*d[j] + d[i]*f[j] + f[i]*d[j])                                                 +
     12.0 *(e[i]*e[j] + f[i]*f[j])                                                                         );
}


void estiva_nsM(MX **Mp, double *S, nde *N)
{
  MX *M; long  a, b, c, A, B, C, e, m, n, I, J, i, j;

  m = dimp2(N); n = dim1(N); initmx(*Mp,m+1,20); M=*Mp;

  for(e=1; e<=n; e++){
    a=N[e].a, b=N[e].b, c=N[e].c, A=N[e].A, B=N[e].B, C=N[e].C;
    setBCD(0, 0, 0, 0, S[e]);
    i=1; 
    foreach(I)&a,&b,&c,&A,&B,&C,end {
      j=1;
      foreach(J)&a,&b,&c,&A,&B,&C,end {
	mx(M,I,J) += l(i,j);
	j++;
      }
      i++;
    }
  }
}
