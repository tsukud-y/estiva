#include "estivaplus.h"

Matrix nsA(Matrix &M,Matrix &K,Matrix &Ax,Matrix &Ay,Matrix &Hx,Matrix  &Hy)
{
  double dt, Re;
  long m = M.capacity();
  Matrix A(m*2+getMatrixsize(Hx));

  dt = tau();
  forMatrixNonzero(M){ 
    A[i][j]     = M[i][j]/dt;
    A[m+i][m+j] = M[i][j]/dt;
  }

  opf(Re,250);
  forMatrixNonzero(K) {
    A[i][j]     += K[i][j]/Re;
    A[m+i][m+j] += K[i][j]/Re;
  }

  forMatrixNonzero(Ax) {
    A[i][j]     += Ax[i][j];
    A[m+i][m+j] += Ax[i][j];
  }

  forMatrixNonzero(Ay) {
    A[i][j]     += Ay[i][j];
    A[m+i][m+j] += Ay[i][j];
  }

  forMatrixNonzero(Hx) {
    A[i][2*m+j] = -Hx[i][j];
    A[2*m+j][i] = -Hx[i][j];
  }

  forMatrixNonzero(Hy) {
    A[m+i][2*m+j] = -Hy[i][j];
    A[2*m+j][m+i] = -Hy[i][j];
  }
  
  return A;
}
