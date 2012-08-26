#include "estivaplus/Mesh.h"

double Mesh::Xmax(vector<Xyc> &Z)
{
  unsigned long i;
  double xmax;

  xmax = Z[1].x;
  for(i=2;i<Z.size()-3;i++)
    if ( xmax < Z[i].x ) xmax = Z[i].x;
  return xmax;
}
