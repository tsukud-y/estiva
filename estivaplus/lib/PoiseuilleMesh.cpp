using namespace std;
#include <estivaplus.h>

void PoiseuilleMesh(double h, double W, double H)
{
  double x, y;

  for ( x = 0.0; x <= W; x+=h) {
    Putnode(x,H,  "G2");
    Putnode(x,0.0,"G2");
  }
  for ( y = 0.0; y <= H; y+=h) {
    Putnode(0.0,y,"G1");
    Putnode(W,  y,"G3");
  }

  for ( x = 0.0; x <= W; x+= h)
    for ( y = 0.0; y <= H; y+= h)
      Putnode(x,y,"");

  Putnode(h/2,h/2,"");
  Putnode(W-h/2,H-h/2,"");

  GenMesh();
}
