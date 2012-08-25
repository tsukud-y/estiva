#include "Mesh.h"

void printZ(vector<Xyc> &Z)
{
  unsigned long i;

  for ( i = 0; i < Z.size(); i++)
    printf("%f %f %s\n",Z[i].x,Z[i].y,Z[i].label.c_str());
}
