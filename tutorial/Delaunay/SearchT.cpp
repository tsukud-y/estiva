#include "Delaunay.h"

static double fplane(double x, double y, double z,
                     double x0,double y0,double z0,
                     double x1,double y1,double z1,
                     double x2,double y2,double z2)
{ double xn,yn,zn;

  xn = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
  yn = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
  zn = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);

  return xn*x+yn*y+zn*z-(xn*x0+yn*y0+zn*z0);
}


long SearchT(vector<Xyc>&Z, vector<Nde>&N, long i)
{
  {
    int e;
    double x, y;
    x = Z[i].x, y = Z[i].y;


    for(e=1;e<N.size();e++){
      long a,b,c; double x0,y0,x1,y1,x2,y2;
      a = N[e].a; b = N[e].b; c = N[e].c;
      x0 = Z[a].x; x1 = Z[b].x; x2 = Z[c].x;
      y0 = Z[a].y; y1 = Z[b].y; y2 = Z[c].y;

      if     (fplane(x,y,0.0, x0,y0,1.0, x1,y1,0.0, x2,y2,0.0)>0.0) continue;
      else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,1.0, x2,y2,0.0)>0.0) continue;
      else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,0.0, x2,y2,1.0)>0.0) continue;
      else return e;
    }
    abort();
  }
}

