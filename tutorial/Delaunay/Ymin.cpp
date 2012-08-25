#include "Mesh.h"

double Mesh::Ymin(vector<Xyc> &Z)
{
  unsigned long i;
  double ymin;

  ymin = Z[1].x;
  for(i=2;i<Z.size()-3;i++)
    if ( ymin > Z[i].x ) ymin = Z[i].x;
  return ymin;
}
