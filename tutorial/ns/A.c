#include "ns.h"


void estiva_nsA(MX **Ap, double *x, double *b, MX *K, MX *M, MX *Hx, MX *Hy, MX *Ax, MX *Ay, double tau, long w)
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
    mx(A,   i,   j) = mx(M,i,j);
    mx(A, m+i, m+j) = mx(M,i,j);
  }

  fornonzeromx(K,i,j) {
    mx(A,   i,   j) += tau*mx(K,i,j);
    mx(A, m+i, m+j) += tau*mx(K,i,j);
  }

  fornonzeromx(Ax,i,j) {
    mx(A,   i,   j) += tau*mx(Ax,i,j);
  }

  fornonzeromx(Ay,i,j) {
    mx(A, m+i, m+j) += tau*mx(Ay,i,j);
  }
  
  fornonzeromx(Hx,i,j) {
    mx(A,    i,2*m+j) = -tau*mx(Hx,i,j);
    mx(A,2*m+j,    i) = -tau*mx(Hx,i,j);
  }

  fornonzeromx(Hy,i,j) {
      mx(A,  m+i,2*m+j) = -tau*mx(Hy,i,j);
      mx(A,2*m+j,  m+i) = -tau*mx(Hy,i,j);
    }
}
