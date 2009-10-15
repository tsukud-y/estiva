#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

void estiva_mulmx(double **tp, MX *A, double *x){
  long i, j, J, m = A->m, n = A->n;
  double *t;
  
  ary1(*tp, m+1);
  
  t  = *tp;
  
  for(i=0; i< m; i++) t[i] = 0.0;

  mx(A,1,1) = mx(A,1,1);

  for(i=0; i< m; i++) for(j=0; j< n; j++) {
      J = A->IA[i][j];
      if (J != 0) t[i] += A->A[i][j]*x[J-1];
  }
}
