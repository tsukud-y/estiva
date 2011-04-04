#include <stdio.h>
#include <string.h>
#include <estiva/op.h>
#include <estiva/solver.h>

int estiva_bicgsolver(void *A, double *x, double *b);
int estiva_gausssolver(void *A, double *x, double *b);
int estiva_blusolver(void *A, double *x, double *b);
int estiva_cgssolver(void *A, double *x, double *b);
int estiva_bicgstabsolver(void *A, double *x, double *b);
int estiva_pcgssolver(void *A, double *x, double *b);    
int estiva_qmrsolver(void *A, double *x, double *b);

int estiva_solver(void *A, double *x, double *b)
{
  if (!strcmp(getop("-solver"),"pcgs") )  return estiva_pcgssolver(A,x,b);
  if (!strcmp(getop("-solver"),"bicg") )  return estiva_bicgsolver(A,x,b);
  if (!strcmp(getop("-solver"),"gauss") ) return estiva_gausssolver(A,x,b);
  if (!strcmp(getop("-solver"),"blu") )   return estiva_blusolver(A,x,b);
  if (!strcmp(getop("-solver"),"cgs") )   return estiva_cgssolver(A,x,b);
  if (!strcmp(getop("-solver"),"bicgstab"))return estiva_bicgstabsolver(A,x,b);
  if (!strcmp(getop("-solver"),"qmr"))     return estiva_qmrsolver(A,x,b);
  return estiva_qmrsolver(A,x,b);
}
