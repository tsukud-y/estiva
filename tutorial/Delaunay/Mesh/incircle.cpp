#include "estivaplus/Mesh.h"

static double incircleVal(int p,int e2,vector<Xyc>&Z,vector<Nde>&N)
{ double a, b, c, x,y,x0,y0,x1,y1,x2,y2;

  x=Z[p].x; y=Z[p].y;
  x0=Z[N[e2].a].x;  y0=Z[N[e2].a].y;
  x1=Z[N[e2].b].x;  y1=Z[N[e2].b].y;
  x2=Z[N[e2].c].x;  y2=Z[N[e2].c].y;

  cramer3(&a,&b,&c, 
          x0,y0,1.0,
          x1,y1,1.0,
          x2,y2,1.0, 
          x0*x0 + y0*y0,
          x1*x1 + y1*y1, 
          x2*x2 + y2*y2);
  return distance2(x0,y0,a/2.0,b/2.0)-distance2(x,y,a/2.0,b/2.0);
}

int Mesh::incircle(int p,int e2,vector<Xyc>&Z,vector<Nde>&N)
{
  return incircleVal(p, e2, Z, N)>0?1:0;
}

int Mesh::incircleDelta(int p,int e2,vector<Xyc>&Z,vector<Nde>&N)
{
  return incircleVal(p, e2, Z, N)>=0?1:0;
}
