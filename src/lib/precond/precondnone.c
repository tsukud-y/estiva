#include <stdio.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>
#include <estiva/op.h>

void estiva_precondnone(long n, double *x, double *b)
{
  long i;
  for (i=0; i<n; i++) x[i] = b[i];
}
