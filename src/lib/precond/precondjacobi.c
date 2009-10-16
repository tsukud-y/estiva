#include <stdio.h>
#include <string.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/precond.h>
#include <estiva/solver.h>
#include <estiva/op.h>

void estiva_precondjacobi(MX *A, double *x, double *D, double *b)
{
  int i,j,J,k,n;
  static double *xold;

  n = dim1(D);

  ary1(xold,n+1);
    
  for (i=0; i<n; i++) xold[i] = 0.0;

  for(k=0;k<3;k++) {
    for(i=0;i<n;i++) {
      x[i]=b[i];
      for(j=0; j< A->n; j++) {
	J = A->IA[i][j];
	if (J != 0) if ( J-1 != i) if(A->A[i][j] !=0.0) {
	  x[i] -= A->A[i][j]*xold[J-1];
	}
      }
      x[i]=x[i]*D[i];
    }
    for(i=0;i<n;i++) xold[i]=x[i];
  }
}

