#include "estivaplus.h"

void nsRhs(Vector &b, Matrix &M, Vector &x)
{
  long m = M.capacity();

  for (long i=0;i<2*m;i++) b[i] = 0.0;

  forMatrixNonzero(M) {
    b[i]   += M[i][j]*x[j];
    b[i+m] += M[i][j]*x[j+m];
  }

  double dt = tau();
  for(long i=0;i<2*m;i++) b[i] /= dt;
}
