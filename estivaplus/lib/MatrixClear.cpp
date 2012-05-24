#include "estivaplus.h"

void MatrixClear(Matrix &A)
{
  unsigned long n = A.capacity();
  
  for (unsigned long i=0; i<n; i++)
    A[i].clear();
  A.clear();
}
