#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mesh.h>

xyc *G_(xyc *Z, nde *N)
{
  static xyc *G;
  double ax,ay,bx,by,cx,cy;
  long i, n;
  n = dim1(N);
  ary1(G,n+1);
  
  for(i=1;i<=n;i++){
    ax = Z[N[i].a].x;
    bx = Z[N[i].b].x;
    cx = Z[N[i].c].x;
    ay = Z[N[i].a].y;
    by = Z[N[i].b].y;
    cy = Z[N[i].c].y;
    G[i].x = (ax+bx+cx)/3.0;
    G[i].y = (ay+by+cy)/3.0;
  }
  return G;
}
