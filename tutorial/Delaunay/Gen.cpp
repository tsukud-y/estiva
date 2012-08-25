#include "Mesh.h"
#undef push
#undef pop
#include <stack>
#include <cmath>
#include <unistd.h>
#include <algorithm>

void Mesh::Gen(vector<Xyc> &Z, vector<Nde> &N)
{
  Mesh::GenSuperNodes(Z,N);
  Mesh::DelaunayAlgo(Z,N);
  Mesh::VanishSuperNodes(Z,N);
  Mesh::VanishBT(Z,N);
  Mesh::SortTri(Z,N);
  Mesh::GenRelation(Z,N);
  Mesh::Normalization(Z,N);
  Mesh::Polynomial2(Z,N);
}
