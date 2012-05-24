#include "estivaplus.h"
#include "stwart.h"

Vector genx(Vector &x0, vector<long>&JCOL)
{
  Vector x(x0.size());

  forVector(x0,i) x[i] = x0[JCOL[i]-1];

  return x;
}
