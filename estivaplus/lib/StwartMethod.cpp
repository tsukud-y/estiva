#include "estivaplus.h"
#include "stwart.h"

void StwartMethod(Matrix &A,vector<long>&JROW,vector<long>&JCOL)
{

  Matrix       AT = TransMatrix(A);

  vector<long> MA = genMA(AT);

  vector<long> IA = genIA(AT,MA);

  MatrixClear(AT);

  Stwart(IA,MA,JROW,JCOL);

  MA.clear();
  IA.clear();
}
