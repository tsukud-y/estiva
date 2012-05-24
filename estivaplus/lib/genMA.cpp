#include "estivaplus.h"
#include "stwart.h"

vector<long> genMA(Matrix &AT)
{
  vector<long> MA(AT.size()+1);

  forMatrix(AT,i,j) MA[i+1]++; 

  MA[0] = 0;
  forVector(MA,i) if(i<MA.size()-1) MA[i+1] += MA[i];

  return MA;
}
