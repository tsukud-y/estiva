#include "estivaplus.h"

void PlotP2(Vector &xv)
{
  static xyc *Z; static nde *N; static double *S;
  static double *x;
  estiva_getZNS(&Z,&N,&S);

  ary1(x,xv.capacity()+1);
  Vector2ary(xv,x);
  estiva_pltp2(x,Z,N);
}
