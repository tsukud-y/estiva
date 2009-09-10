#include <stdio.h>
#include <estiva/op.h>
#include <estiva/solver.h>

int estiva_solver(void *A, double *x, double *b)
{

  if (!strcmp(getop("-solver"),"gauss") ) return estiva_gausssolver(A,x,b);
  if (!strcmp(getop("-solver"),"blu") ) return estiva_blusolver(A,x,b);
  return estiva_pcgssolver(A,x,b);
}
