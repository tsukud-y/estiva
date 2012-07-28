#include "Delaunay.h"

void Putmesh(vector<Xyc> &Z, vector<Nde> &N)
{
  printf("<xyc>\n");
  for (unsigned long i=0; i<Z.size(); i++)
    printf("%f %f %s\n",Z[i].x,Z[i].y,Z[i].label.c_str());

  printf("<nde>\n");
  for (unsigned long i=0; i<N.size(); i++)
    printf("%ld %ld %ld %ld %ld %ld\n",
           N[i].a,N[i].b,N[i].c,N[i].A,N[i].B,N[i].C);

 
}
