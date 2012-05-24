#include "estivaplus.h"
#include "stwart.h"

Vector makeb(long M)
{
  Vector b(M);

  forVector(b,i)  b[i] = (double)(i+1)*(i+1);

  return b;
}
