#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/que.h"
#include "estiva/op.h"


void estiva_genmesh(void *p, xyc **Zp, nde **Np)
{

  static xyc *Z, *ep, e;
  static nde *N;
  long        i,  n;
  que **qp;
  qp = p;

  n = 0; forq(*qp,e) n++;
  
  ary1(Z,n+4);
  i = n; forq(*qp,ep) {
    Z[i].x = ep->x, Z[i].y =ep->y, Z[i].label = ep->label;
    i--;
  }
  while (pop(*qp,e));

  if (defop("-xmesh")) xmesh(Z);
  delaunay(Z,N);
  *Zp = Z, *Np = N;
}
