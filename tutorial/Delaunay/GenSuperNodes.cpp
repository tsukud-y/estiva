#include "Mesh.h"

void Mesh::GenSuperNodes(vector<Xyc>&Z, vector<Nde>&N) {  
  Xyc z;
  z.x=0.0, z.y=0.0, z.label="super_node"; Z[0] = z;
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);
  z.x=0.0, z.y=0.0, z.label="super_node"; Z.push_back(z);

  unsigned long n = Z.size();

  Z[n-3].x=Mesh::Xmin(Z), Z[n-3].y=Mesh::Ymax(Z);
  Z[n-2].x=Mesh::Xmax(Z), Z[n-2].y=Mesh::Ymax(Z);
  Z[  0].x=Mesh::Xmin(Z), Z[  0].y=Mesh::Ymin(Z);
  Z[n-1].x=Mesh::Xmax(Z), Z[n-1].y=Mesh::Ymin(Z);

  double length 
    =sqrt(pow(Mesh::Xmax(Z)-Mesh::Xmin(Z),2)+pow(Mesh::Ymax(Z)-Mesh::Ymin(Z),2))/2;

  Z[n-3].x-=length, Z[n-3].y+=length;
  Z[n-2].x+=length, Z[n-2].y+=length;
  Z[n-1].x+=length, Z[n-1].y-=length;
  Z[  0].x-=length, Z[  0].y-=length;

  Nde nd;
 
  nd.a=  0, nd.b=0,   nd.c=0,   nd.A=0, nd.B=0, nd.C=0;
  N.push_back(nd);
  nd.a=  0, nd.b=n-1, nd.c=n-3, nd.A=2, nd.B=0, nd.C=0;
  N.push_back(nd);
  nd.a=n-2, nd.b=n-3, nd.c=n-1, nd.A=1, nd.B=0, nd.C=0;
  N.push_back(nd);
}
