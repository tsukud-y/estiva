#include "estivaplus.h"
#include "stwart.h"

vector<long> genIA(Matrix &AT, vector<long>&MA)
{
  vector<long> IA(MA[AT.size()]);

  long k = 0;
  forMatrix(AT,i,j) IA[k++] = j+1;

  return IA;
}
