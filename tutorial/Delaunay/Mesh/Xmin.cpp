#include "estivaplus/Mesh.h"

double Mesh::Xmin(vector<Xyc> &Z)
{
  unsigned long i;
  double xmin;

  xmin = Z[1].x;
  for(i=2;i<Z.size()-3;i++)
    if ( xmin > Z[i].x ) xmin = Z[i].x;
  return xmin;
}
