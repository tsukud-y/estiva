#include "fem.h"

#include "estiva/ary.h"


void estiva_femdelta(double **Sp, xyc *Z, nde *N)
{
  static double *S;
  long i, e;
  double xi, xj, xk, yi, yj, yk;

  e = dim1(N);
  ary1(S,e+1);

  for( i=1; i<=e; i++ ){
    xi = Z[N[i].a].x;
    xj = Z[N[i].b].x;
    xk = Z[N[i].c].x;
    yi = Z[N[i].a].y;
    yj = Z[N[i].b].y;
    yk = Z[N[i].c].y;
    S[i] = xi*yj+xj*yk+xk*yi-yi*xj-yj*xk-yk*xi;
    S[i] /= 2.0;
  }
  *Sp = S;
}
