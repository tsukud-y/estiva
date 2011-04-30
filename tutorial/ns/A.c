#include "ns.h"


void estiva_nsA(MX **Ap, double *x, double *b, MX *K, MX *M, MX *Hx, MX *Hy, MX *Ax, MX *Ay, long w)
{
  static MX *A;
  long   i, j, NUM, m, n;
  static xyc *Z; static nde *N; static double *S;

  getZNS(Z,N,S);
  m   = M->n;
  n   = dim1(Z);
  NUM = m*2+n;
  initmx(*Ap, NUM+1, w); A = *Ap;
  clearmx(A);

  fornonzeromx(M,i,j) {
    mx(A,   i,   j) = mx(M,i,j)/tau();
    mx(A, m+i, m+j) = mx(M,i,j)/tau();
  }

  fornonzeromx(K,i,j) {
    mx(A,   i,   j) += mx(K,i,j)/Re();
    mx(A, m+i, m+j) += mx(K,i,j)/Re();
  }

  fornonzeromx(Ax,i,j) {
    mx(A,   i,   j) += mx(Ax,i,j);
    mx(A, m+i, m+j) += mx(Ax,i,j);
  }

  fornonzeromx(Ay,i,j) {
    mx(A,   i,   j) += mx(Ay,i,j);
    mx(A, m+i, m+j) += mx(Ay,i,j);
  }
  
  fornonzeromx(Hx,i,j) {
    mx(A,    i,2*m+j) = -mx(Hx,i,j);
    mx(A,2*m+j,    i) = -mx(Hx,i,j);
  }

  fornonzeromx(Hy,i,j) {
      mx(A,  m+i,2*m+j) = -mx(Hy,i,j);
      mx(A,2*m+j,  m+i) = -mx(Hy,i,j);
    }
}
