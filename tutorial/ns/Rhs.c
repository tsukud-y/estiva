#include "ns.h"


void estiva_nsRhs(double *b, MX *M, double *x)
{
  long   i, m;
  static double *bx, *by;

  m = M->n;

  mulmx(bx,M,x+1);
  mulmx(by,M,x+m+1);

  for(i=1;i<=m;i++){
    b[i] = bx[i-1];
    b[i+m] = by[i-1];
  }
}
