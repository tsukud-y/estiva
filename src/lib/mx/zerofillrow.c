#include "estiva/mx.h"


void estiva_zerofillrow(MX *A, long i)
{
  long j;
  mx(A,1,1) = mx(A,1,1);

  i--;
  for ( j = 0; j < A->w; j++ )
    if(A->A[i][j] != 0.0)
      A->A[i][j] = 0.0;
}
