#include <stdio.h>
#include <estiva/op.h>
#include "msh.h"
#include "ary.h"
#include "spm.h"

#define A(i,j) mx(A,i,j)

static void cplot(xyc *Z,nde *N,void* A)
{
  int i,j,n;
  n = dim1(Z);
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++) 
      if(A(i,j) == 1.0) printf("£±");
      else if(A(i,j) != 0.0) printf("¢£");
      else  printf("¢¢");
    printf("\n");
  }
}
static void splot(xyc *Z, nde *N,double* u)
{
  int e,a,b,c;
  for(e=1;e<=dim1(N);e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    printf("%f %f %f\n",Z[a].y,Z[a].x,-u[a]);
    printf("%f %f %f\n",Z[b].y,Z[b].x,-u[b]);
    printf("%f %f %f\n",Z[c].y,Z[c].x,-u[c]);
    printf("%f %f %f\n",Z[a].y,Z[a].x,-u[a]);
    printf("\n\n");
  }
}
static void makeuv(xyc *Z, nde *N,double *S, double *b, double *u, double *v)
{
  long e, n, i, j, k;
  double x1, x2, x3, y1, y2, y3;
  n = dim1(N);
  for(e=1; e<=n; e++){
    i = N[e].a,  j = N[e].b,  k = N[e].c;
    x1 = Z[i].x,  y1 = Z[i].y;
    x2 = Z[j].x,  y2 = Z[j].y;
    x3 = Z[k].x,  y3 = Z[k].y;
    u[e] =  (y1*b[i] + y2*b[j] + y3*b[k])/2.0/S[e];
    v[e] = -(x1*b[i] + x2*b[j] + x3*b[k])/2.0/S[e];
  }
}
static void uvplot(xyc *Z, nde *N, double *u, double *v)
{
  static double *Gx, *Gy;
  long e, n, a, b, c;
  double xa, xb, xc, ya, yb, yc, t;
  n = dim1(N);
  ary1(Gx,n+1),  ary1(Gy,n+1);
  for(e=1; e<=n; e++){
    a = N[e].a,  b = N[e].b,  c = N[e].c;
    xa = Z[a].x,      xb = Z[b].x,      xc = Z[c].x;
    ya = Z[a].y,      yb = Z[b].y,      yc = Z[c].y;
    Gx[e] = (xa + xb + xc)/3.0;
    Gy[e] = (ya + yb + yc)/3.0;
  }
  t = 0.03;
  for(e=1; e<=n; e++){
    u[e] *= 0.1;
    v[e] *= 0.1;
    printf("%f %f\n",Gx[e],Gy[e]);
    printf("%f %f\n",Gx[e]+t*u[e],Gy[e]+t*v[e]);
    printf("\n");
  }
}
void plot(xyc *Z, nde *N, long *A, double *b)
{
  static double *S;
  S = S_(Z,N);
  if(defop("-uv")){
    static double *u, *v;
    ary1(u,dim1(N)+1);
    ary1(v,dim1(N)+1);
    makeuv(Z,N,S,b,u,v);
    uvplot(Z,N,u,v);
  }
  else if(defop("-cplot")) cplot(Z,N,A);
  else splot(Z, N, b);
}

