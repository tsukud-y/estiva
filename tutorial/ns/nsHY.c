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


static double hy(long i, long j)
{
  return (Delta/12.0) * (12.0*(alphaC(i)*ad[j])                                                       
			 +4.0*( betaC(i)*ad[j] + alphaC(i)*bd[j] + gammaC(i)*ad[j] + alphaC(i)*cd[j]) 
			 +1.0*( betaC(i)*cd[j] + gammaC(i)*bd[j])                                     
			 +2.0*( betaC(i)*bd[j] + gammaC(i)*cd[j])                                     );
}


void estiva_HY(MX **Hyp, double *S, xyc *Z, nde *N)
{
  MX *Hy; long a, b, c, A, B, C, e, m, n, I, J, i, j; double s;
  m = dimp2(N); n = dim1(N); initmx(*Hyp,m+1,21); Hy = *Hyp;
  
  for(e=1;e<=n;e++){
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C; s=S[e];
    setBCD(NAN, NAN, (Z[c].x-Z[b].x)/(2.0*s), (Z[a].x-Z[c].x)/(2.0*s), s); 
    i=1; 
    foreach(I)&a,&b,&c,&A,&B,&C,end {
      j=1;
      foreach(J)&a,&b,&c,end {
        mx(Hy,I,J) += hy(i,j);
        j++;
      }
      i++;
    }
  }
}
