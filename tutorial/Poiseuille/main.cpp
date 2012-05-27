#include "estivaplus.h"

void PoiseuilleMesh()
{
  double x, y, h=0.1, W=6, H=1;

  for ( x = h; x <= W-h; x+= h)
    for ( y = h; y <= H-h; y+= h)
      putnode(x,y,NULL);

  for ( x = 0.0; x < W+h; x+=h) {
    putnode(x,H,  "G2");
    putnode(x,0.0,"G2");
  }
  for ( y = 0.0; y < H+h; y+=h) {
    putnode(0.0,y,"G1");
    putnode(W,  y,"G3");
  }
  GenerateMesh();
}

int main(int argc, char **argv){
  initop(argc,argv);
  PoiseuilleMesh();
  return 0;
}
