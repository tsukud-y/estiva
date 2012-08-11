#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>

int main(int argc, char ** argv)
{
  vector<Xyc> Z; vector<Nde> N;  
  Xyc z;

  Z.push_back(z);

  for ( z.y = 0.0; z.y <= 1.0; z.y+=0.125)
    for ( z.x = 0.0; z.x <= 1.0; z.x+=0.125)
      {
	if ( z.x == 0.0 || z.y == 0.0 || z.x == 1.0 || z.y == 1.0)
	  z.label="G";
	else z.label="";
	Z.push_back(z);
      }

  FILE *pp = popen("gnuplot","w");

  GenSuperNodes(Z,N);
  DelaunayAlgo(Z,N);
  VanishSuperNodes(Z,N);
  VanishBT(Z,N);
  SortTri(Z,N);
  GenRelation(Z,N);
  Normalization(Z,N);

  Putmesh(Z,N);
  Xmesh(pp,Z,N);
  sleep(300);
}
