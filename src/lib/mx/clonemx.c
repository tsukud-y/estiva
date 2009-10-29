#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

void estiva_clonemx(MX **MT, MX *M)
{
  long i, j, J;  

  initmx((*MT), M->m+1, M->n+1);

  mx(M,1,1) = mx(M,1,1);

  for(i=1;i<= M->m;i++) for(j=0; j< M->n; j++) {
      J = M->IA[i-1][j];
      if (J != 0) mx((*MT),i,J) = M->A[i-1][j];
  }
}

