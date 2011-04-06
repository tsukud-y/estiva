#include <stdio.h>
#include <estiva/ary.h>
#include <estiva/mx.h>
#include <estiva/solver.h>

void estiva_transmx(MX **MT, MX *M)
{
  long i, j, J;  

  initmx((*MT), M->n+1, M->w+1);

  mx(M,1,1) = mx(M,1,1);

  for(i=1;i<= M->n;i++) for(j=0; j< M->w; j++) {
      J = M->IA[i-1][j];
      if (J != 0) mx((*MT),J,i) = M->A[i-1][j];
  }
}

