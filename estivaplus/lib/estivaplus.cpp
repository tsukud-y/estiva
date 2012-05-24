#include "estivaplus.h"

std::map<unsigned int, double>::iterator estiva_aiterator;

Vector ary2Vector(double *x)
{
  long i,NUM;
  NUM = dim1(x);
  Vector xv(NUM);
  for (i=0; i<NUM; i++) xv[i] = x[i+1];
  return xv;
}

void Vector2ary(Vector &bv, double *b)
{
  long i, NUM;
  NUM = dim1(b);
  for (i=0; i<NUM; i++) b[i+1] = bv[i];
}

Matrix mx2Matrix(MX *A)
{
  long i,j,n;
  n = A->n;

  Matrix a(n);
  a.clear();

  fornonzeromx(A,i,j) {
    a[i-1][j-1] = (*estiva_mx(A,i,j));
  }

  return a;
}

void Matrix2mx(Matrix &a, MX **Ap)
{
  static MX  *A;
  unsigned int i,j,maxj = 0;
  for (i=0; i<a.capacity(); i++) {
    std::map<unsigned int, double>::iterator ai = a[i].begin();
    unsigned int j = 0;
    while ( ai != a[i].end()){
      ++j;
      ++ai;
    }
    if ( j > maxj) maxj = j;
  }
  estiva_initmx(Ap, a.capacity()+1, maxj+1); A = *Ap;
  estiva_clearmx(A);

  for ( i=0; i<a.capacity(); i++) {
    std::map<unsigned int, double>::iterator ai = a[i].begin();
    while ( ai != a[i].end()) {
      j = (*ai).first;
      (*estiva_mx(A,i+1,j+1)) = a[i][j];
      ++ai;
    }
  }
}


long getMatrixsize(Matrix &Hxm)
{
  unsigned int i, j, jmax;
  for ( jmax=0, i=0; i<Hxm.capacity(); i++ ) {
    std::map<unsigned int, double>::iterator ai = Hxm[i].begin();
    while ( ai != Hxm[i].end()) {
      j = (*ai).first;
      if ( jmax < j ) jmax = j;
      ++ai;
    }
  }
  return jmax+1;
}


