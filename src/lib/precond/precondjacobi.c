#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/precond.h>
#include <estiva/solver.h>
#include <estiva/op.h>

void estiva_precondjacobi(MX *A, double *x, double *D, double *b)
{
  int i,j,J,k,n, step = 2;
  static double *xo;

  if ( defop("-jacobistep") ) step = atoi(getop("-jacobistep"));

  n = dim1(D);

  ary1(xo,n+1);
    
  for (i=0; i<n; i++) xo[i] = 0.0;

  for(k=0;k<step;k++) {
    for(i=0;i<n;i++) {
      x[i]=b[i];
      for(j=0; j< A->w; j++) {
	J = A->IA[i][j];
	if (J != 0) if ( J-1 != i) if(A->A[i][j] !=0.0) {
	  x[i] -= A->A[i][j]*xo[J-1];
	}
      }
      x[i]=x[i]*D[i];
    }
    for(i=0;i<n;i++) xo[i]=x[i];
  }
}
