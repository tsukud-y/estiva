#include <stdio.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/esolver.h>

extern double estiva_maxesolver(double **A, double *x);
extern double estiva_minesolver(double **A, double *x);

double estiva_esolver(double **A, double *x)
{
  if (!strcmp(getop("-esolver"),"min"))
    return estiva_minesolver(A,x);
  return estiva_maxesolver(A,x);
}
