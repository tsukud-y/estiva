#include "ns.h"
#include "fem.h"


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/op.h"
#include "estiva/mx.h"
#include "estiva/foreach.h"
#include "estiva/que.h"
#include "estiva/solver.h"


void estiva_nsA(MX **Ap, double *x, double *b, xyc *Z, nde *N, MX *K, MX *M, MX *Hx, MX *Hy, MX *AX, MX *AY, double tau)
{
  static MX *A;
  long   i, j, NUM, m, n;

  m   = dimp2(N);
  n   = dim1(Z);
  NUM = m*2+n;
  initmx(*Ap, NUM+1, 50); A = *Ap;

  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= m; j++ ) {
      mx(A,  i,   j) = mx(M,i,j) + tau*mx(K,i,j) + tau*mx(AX,i,j);
      mx(A,m+i, m+j) = mx(M,i,j) + tau*mx(K,i,j) + tau*mx(AY,i,j);
    }
  for ( i = 1; i <= m; i++ ) for ( j = 1; j <= n; j++ ) {
      mx(A,    i,2*m+j) = -tau*mx(Hx,i,j);
      mx(A,2*m+j,    i) = -tau*mx(Hx,i,j);
      mx(A,  m+i,2*m+j) = -tau*mx(Hy,i,j);
      mx(A,2*m+j,  m+i) = -tau*mx(Hy,i,j);
    }
}
