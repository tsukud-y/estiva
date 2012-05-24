#include "estivaplus.h"
#include "stwart.h"

void SolverWithStwart(Matrix &A, Vector &x, Vector &b)
{
  long   M  = A.size();
  vector<long> JROW(M), JCOL(M);
  StwartMethod(A,JROW,JCOL);
  Matrix A0 = genA0(A,JROW,JCOL);
  Vector b0 = genb0(b,JROW);

  Solverorg(A0,x,b0);

  x = genx(x,JCOL);
}
