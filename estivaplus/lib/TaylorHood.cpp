#include "estivaplus.h"

Matrix TaylorHood_M()
{
  static MX *M;
  estiva_TaylorHood_M(&M,36);
  Matrix Mm = mx2Matrix(M);
  return Mm;
}

Matrix TaylorHood_K()
{
  static MX *K;
  estiva_TaylorHood_K(&K,32);
  Matrix Km = mx2Matrix(K);
  return Km;
}

Matrix TaylorHood_Hx()
{
  static MX *Hx;
  estiva_TaylorHood_Hx(&Hx,15);
  Matrix Hxm = mx2Matrix(Hx);
  return Hxm;
}

Matrix TaylorHood_Hy()
{
  static MX *Hy;
  estiva_TaylorHood_Hy(&Hy,15);
  Matrix Hym = mx2Matrix(Hy);
  return Hym;
}

Matrix TaylorHood_Ax(Vector &xv)
{
  static MX *Ax;
  static double *x;
  ary1(x,xv.capacity()+1);
  Vector2ary(xv,x);
  estiva_TaylorHood_Ax(&Ax,x,34);
  Matrix Axm = mx2Matrix(Ax);
  return Axm;
}

Matrix TaylorHood_Ay(Vector &xv)
{
  static MX *Ay, *M;
  static double *x;
  ary1(x,xv.capacity()+1);
  Vector2ary(xv,x);
  estiva_TaylorHood_M(&M,36);
  estiva_TaylorHood_Ay(&Ay,x+M->n,34);
  Matrix Aym = mx2Matrix(Ay);
  return Aym;
}
