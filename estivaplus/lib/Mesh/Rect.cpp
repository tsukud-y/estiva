#include "estivaplus/Mesh.h"

void Mesh::Rect(vector<Xyc>&Z, vector<Nde>&N)
{
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

  Mesh::Gen(Z,N);
}
