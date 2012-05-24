#include "estivaplus.h"
#include "stwart.h"

Matrix genA0(Matrix &A,vector<long>&JROW,vector<long>&JCOL)
{
  long M = A.size();
  Matrix B(M), B2(M);
  forMatrix(A,i,j) B[JROW[i]-1][j]  = A[i][j];
  forMatrix(B,i,j) B2[i][JCOL[j]-1] = B[i][j];
  return B2;
}
