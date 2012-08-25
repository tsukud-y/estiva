#include <estivaplus.h>

#include <algorithm>
#include <cmath>

using namespace std;

struct Xyc { double x, y; string label; };

bool LessXY(const Xyc& rLeft, const Xyc& rRight) { 

  if ( rLeft.x != rRight.x)
    return rLeft.x < rRight.x; 

  return rLeft.y < rRight.y;
}

vector<Xyc> Zv;

void Putnode(double x, double y, string label)
{
  Xyc z;
  unsigned long i;
  
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
  unsigned long i;
  for ( i = 0; i < Zv.size(); i++)
    printf("%f %f %s\n",Zv[i].x, Zv[i].y, Zv[i].label.c_str());
}


void PushXyc(void)
{
  unsigned long i;
  for ( i = 0; i < Zv.size(); i++) {
    if (Zv[i].label == "")
      putnode(Zv[i].x, Zv[i].y, NULL);
    else
      putnode(Zv[i].x, Zv[i].y, Zv[i].label.c_str());
  }
}

void GenMeshP(void){
  
  SortMesh();
  //PutXyc();
  PushXyc();
  GenerateMesh();
}

