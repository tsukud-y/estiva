#include "Mesh.h"
#include <cmath>

void Mesh::GenSuperNodes(vector<Xyc>&Z, vector<Nde>&N) {  
  Xyc z;
  z.x=0.0, z.y=0.0, z.label="super_node"; Z[0] = z;
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);

  unsigned long n = Z.size();

  Z[n-3].x=Mesh::Xmin(Z), Z[n-3].y=Mesh::Ymax(Z);  Z[n-2].x=Mesh::Xmax(Z), Z[n-2].y=Mesh::Ymax(Z);
  Z[  0].x=Mesh::Xmin(Z), Z[  0].y=Mesh::Ymin(Z);  Z[n-1].x=Mesh::Xmax(Z), Z[n-1].y=Mesh::Ymin(Z);

  double length = sqrt( pow(Mesh::Xmax(Z)-Mesh::Xmin(Z),2) + pow(Mesh::Ymax(Z)-Mesh::Ymin(Z),2))/2;

  Z[n-3].x-=length, Z[n-3].y+=length;  Z[n-2].x+=length, Z[n-2].y+=length;
  Z[  0].x-=length, Z[  0].y-=length;  Z[n-1].x+=length, Z[n-1].y-=length;

  Nde nde;
 
  nde.a=  0, nde.b=0,   nde.c=0,   nde.A=0, nde.B=0, nde.C=0;
  N.push_back(nde);
  nde.a=  0, nde.b=n-1, nde.c=n-3, nde.A=2, nde.B=0, nde.C=0;
  N.push_back(nde);
  nde.a=n-2, nde.b=n-3, nde.c=n-1, nde.A=1, nde.B=0, nde.C=0;
  N.push_back(nde);
}
