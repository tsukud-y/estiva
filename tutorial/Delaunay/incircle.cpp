#include "Delaunay.h"

#define sarrus(a11,a12,a13,a21,a22,a23,a31,a32,a33) \
  ((a11)*(a22)*(a33)+(a21)*(a32)*(a13)+(a31)*(a12)*(a23)\
   -(a13)*(a22)*(a31)-(a23)*(a32)*(a11)-(a33)*(a12)*(a21))

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
static void cramer3(double *px,double *py,double *pz, 
                    double a11,double a12,double a13, 
                    double a21,double a22,double a23, 
                    double a31,double a32,double a33,
                    double b1, double b2, double b3 )
{ double det;

  *px = sarrus(b1 ,a12,a13, b2 ,a22,a23, b3 ,a32,a33);
  *py = sarrus(a11,b1 ,a13, a21,b2 ,a23, a31,b3 ,a33);
  *pz = sarrus(a11,a12,b1 , a21,a22,b2 , a31,a32,b3 );
  det = sarrus(a11,a12,a13, a21,a22,a23, a31,a32,a33);
  if(det != 0.0){ *px/=det;*py/=det;*pz/=det;}
}



int incircle(int p,int e2,vector<Xyc>&Z,vector<Nde>&N)
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

  return distance2(x0,y0,a/2.0,b/2.0)-distance2(x,y,a/2.0,b/2.0)>0.0? 1:0;
}
