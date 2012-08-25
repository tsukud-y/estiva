#include "Delaunay.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>

void GenMesh(vector<Xyc> &Z, vector<Nde> &N)
{
  GenSuperNodes(Z,N);
  DelaunayAlgo(Z,N);
  VanishSuperNodes(Z,N);
  VanishBT(Z,N);
  SortTri(Z,N);
  GenRelation(Z,N);
  Normalization(Z,N);
  Polynomial2(Z,N);
}
