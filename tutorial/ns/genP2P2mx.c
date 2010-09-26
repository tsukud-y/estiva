#include "ns.h"
#include "fem.h"


#include "estiva/mx.h"
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"

static xyc *Zg;
static nde *Ng;
static double *Sg;

void setZNS(xyc *Z, nde *N, double *S)
{
  Zg = Z; Ng = N; Sg = S;
}

void genP2P2mx(MX **Mp, double (*func)(long i, long j))
{
  MX *M; long  a, b, c, A, B, C, e, m, n, I, J, i, j; double s;
  xyc *Z; nde *N; double *S;
  Z = Zg; N = Ng; S = Sg;

  m = dimp2(N); n = dim1(N); initmx(*Mp,m+1,28); M = *Mp;

  for (e=1; e<=n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C; s=S[e];
    setBCD((Z[b].y-Z[c].y)/(2.0*s), (Z[c].y-Z[a].y)/(2.0*s), (Z[c].x-Z[b].x)/(2.0*s), (Z[a].x-Z[c].x)/(2.0*s), s);
    i=1; 
    foreach(I)&a,&b,&c,&A,&B,&C,end {
      j=1;
      foreach(J)&a,&b,&c,&A,&B,&C,end {
        mx(M,I,J) += func(i,j);
        j++;
      }
      i++;
    }
  }
}

void genP2P1mx(MX **Mp, double (*func)(long i, long j))
{
  MX *M; long  a, b, c, A, B, C, e, m, n, I, J, i, j; double s;
  xyc *Z; nde *N; double *S;
  Z = Zg; N = Ng; S = Sg;

  m = dimp2(N); n = dim1(N); initmx(*Mp,m+1,28); M = *Mp;

  for (e=1; e<=n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C; s=S[e];
    setBCD((Z[b].y-Z[c].y)/(2.0*s), (Z[c].y-Z[a].y)/(2.0*s), (Z[c].x-Z[b].x)/(2.0*s), (Z[a].x-Z[c].x)/(2.0*s), s);
    i=1; 
    foreach(I)&a,&b,&c,&A,&B,&C,end {
      j=1;
      foreach(J)&a,&b,&c,end {
        mx(M,I,J) += func(i,j);
        j++;
      }
      i++;
    }
  }
}
