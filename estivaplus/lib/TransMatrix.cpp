#include "estivaplus.h"
#include "stwart.h"

Matrix TransMatrix(Matrix &A)
{
  Matrix AT(A.size());

  forMatrix(A,i,j) AT[j][i] = A[i][j];

  return AT;
}
