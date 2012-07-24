#include "estivaplus.h"

#include <algorithm>
#include <cmath>

using namespace std;

struct Xyc { double x, y; string label; };

bool LessX(const Xyc& rLeft, const Xyc& rRight) { return rLeft.x < rRight.x; }
bool LessY(const Xyc& rLeft, const Xyc& rRight) { return rLeft.y < rRight.y; }
bool LessXY(const Xyc& rLeft, const Xyc& rRight) { 

  if ( rLeft.x != rRight.x)
    return rLeft.x < rRight.x; 

  return rLeft.y < rRight.y;
}

vector<Xyc> Zv;

void Putnode(double x, double y, string label)
{
  Xyc z;
  long i;
  
  z.x = x, z.y = y; z.label = label;


  for ( i = 0; i < Zv.size(); i++){
    if ( fabs(Zv[i].x-x) < 0.001 && fabs(Zv[i].y-y) < 0.001)
      return;
  }
  Zv.push_back(z);
}

void SortMesh(void)
{
  sort(Zv.begin(), Zv.end(), LessXY);  
}

void PutXyc(void)
{
  long i;
  for ( i = 0; i < Zv.size(); i++)
    printf("%f %f %s\n",Zv[i].x, Zv[i].y, Zv[i].label.c_str());
}


void PushXyc(void)
{
  long i;
  for ( i = 0; i < Zv.size(); i++) {
    if (Zv[i].label == "")
      putnode(Zv[i].x, Zv[i].y, NULL);
    else
      putnode(Zv[i].x, Zv[i].y, Zv[i].label.c_str());
  }
}

void GenMesh(void){
  
  SortMesh();
  //PutXyc();
  PushXyc();
  GenerateMesh();
}

void PoiseuilleMesh()
{
  double x, y, h=0.2, W=6, H=1;

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

int main(int argc, char **argv){
  initop(argc,argv);
  PoiseuilleMesh();
  return 0;
}
