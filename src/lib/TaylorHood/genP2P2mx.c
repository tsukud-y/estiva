#include <stdio.h>
#include "estiva/mx.h"
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"
#include "estiva/TaylorHood.h"

static xyc *Zg;
static nde *Ng;
static double *Sg;

static void setZNS(xyc *Z, nde *N, double *S)
{
  Zg = Z; Ng = N; Sg = S;
}

void estiva_getZNS(xyc **Zp, nde **Np, double **Sp)
{
  *Zp = Zg; *Np = Ng; *Sp = Sg;
}

double *estiva_setmesh(xyc *Z, nde *N)
{
  static double *S;
  femdelta(S,Z,N);
  setZNS(Z,N,S);
  return S;
}


void estiva_genP2P2mx(MX **Mp, double (*func)(long i, long j), long w)
{
  MX *M; long  a, b, c, A, B, C, e, m, n, I, J, i, j; double s;
  xyc *Z; nde *N; double *S;
  Z = Zg; N = Ng; S = Sg;

  m = dimp2(N); n = dim1(N); initmx(*Mp,m+1,w); M = *Mp;

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

void estiva_genP2P1mx(MX **Mp, double (*func)(long i, long j), long w)
{
  MX *M; long  a, b, c, A, B, C, e, m, n, I, J, i, j; double s;
  xyc *Z; nde *N; double *S;
  Z = Zg; N = Ng; S = Sg;

  m = dimp2(N); n = dim1(N); initmx(*Mp,m+1,w); M = *Mp;

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
