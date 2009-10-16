#include <stdio.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>
#include <estiva/op.h>

void estiva_precondscaling(double *x, double *D, double *b)
{
  long i, n;
  n = dim1(D);
  for (i=0; i<n; i++) x[i] = D[i]*b[i];
}
