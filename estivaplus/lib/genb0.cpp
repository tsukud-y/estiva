#include "estivaplus.h"
#include "stwart.h"

Vector genb0(Vector &b, vector<long> &JROW)
{
  Vector b0(b.size());

  forVector(b,i) b0[JROW[i]-1] = b[i];

  return b0;
}
