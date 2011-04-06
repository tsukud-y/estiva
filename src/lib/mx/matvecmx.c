#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

void estiva_matvecmx(MX *A, double *alpha, double *x, double *beta, double *y)
{
  static double *t;
  long i, n = A->n; 

  mulmx(t,A,x);
  for(i=0;i<n;i++) t[i] *= *alpha;
  for(i=0;i<n;i++) t[i] += *beta*y[i];
  for(i=0;i<n;i++) y[i] = t[i];
}
