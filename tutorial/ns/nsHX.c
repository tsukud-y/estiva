#include "ns.h"
#include "fem.h"

#include <math.h>
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/mx.h"

static double ad[] = { NAN,         0.0,  0.0,  1.0, NAN};
static double bd[] = { NAN,         1.0,  0.0, -1.0, NAN};
static double cd[] = { NAN,         0.0,  1.0, -1.0, NAN};

static double hx(long i, long j)
{
  return (Delta/12.0) * (12.0*(alphaB(i)*ad[j])                                                      
			 +4.0*( betaB(i)*ad[j] + alphaB(i)*bd[j] + gammaB(i)*ad[j] + alphaB(i)*cd[j]) 
			 +1.0*( betaB(i)*cd[j] + gammaB(i)*bd[j])                                    
			 +2.0*( betaB(i)*bd[j] + gammaB(i)*cd[j])                                     );
}


void estiva_HX(MX **Hxp, double *S, xyc *Z, nde *N)
{
  MX *Hx;long a, b, c, A, B, C, e, m, n, I, J, i, j; double s;
  m = dimp2(N); n = dim1(N); initmx(*Hxp,m+1,21); Hx = *Hxp;

  for(e=1;e<=n;e++){
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C; s=S[e];
    setBCD((Z[b].y-Z[c].y)/(2.0*s), (Z[c].y-Z[a].y)/(2.0*s), NAN, NAN, s);
    i=1; 
    foreach(I)&a,&b,&c,&A,&B,&C,end {
      j=1;
      foreach(J)&a,&b,&c,end {
        mx(Hx,I,J) += hx(i,j);
        j++;
      }
      i++;
    }
  }
}
