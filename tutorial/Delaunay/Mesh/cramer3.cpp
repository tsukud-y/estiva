#include "estivaplus/Mesh.h"

void cramer3(double *px,double *py,double *pz, 
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


