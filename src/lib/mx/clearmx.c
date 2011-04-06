#include <stdio.h>
#include "estiva/mx.h"

void estiva_clearmx(MX *A){
  long i,j;
  mx(A,1,1) = 0.0;

  for ( i = 0; i < A->n; i++ )
    for ( j = 0; j < A->w; j++ )
      if(A->A[i][j] != 0.0)
        A->A[i][j] = 0.0;
}
