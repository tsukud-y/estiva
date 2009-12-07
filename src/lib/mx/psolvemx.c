#include <stdio.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/precond.h>
#include <estiva/solver.h>
#include <estiva/op.h>

void estiva_psolvemx(MX *A, CRS *pivot, MX *LU, double *D, double *x, double *b)
{
  if (!strcmp(getop("-precond"),"none"))
    precondnone(A->m,x,b); 
  else if (!strcmp(getop("-precond"),"ILU"))
    precondILU(pivot,LU,x,b);
  else if (!strcmp(getop("-precond"),"scaling"))
    precondscaling(x,D,b); 
  else 
    precondjacobi(A,x,D,b);
}
