#include "estivaplus/Mesh.h"

double Mesh::Ymax(vector<Xyc> &Z)
{
  unsigned long i;
  double Ymax;

  Ymax = Z[1].x;
  for(i=2;i<Z.size()-3;i++)
    if ( Ymax < Z[i].x ) Ymax = Z[i].x;
  return Ymax;
}
