#include <stdio.h>
#include "estiva/std.h"
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/solver.h"

void estiva_slimupmx(MX **M, MX *A)
{
  long i, j, J,  maxw=0, ii=0, w=0;  

  mx(A,1,1) = mx(A,1,1);

  for (i = ii; i < A->m; i++, ii++, w=0) 
    for (j = 0; j< A->n; j++) 
      if (A->IA[i][j] != 0.0) {
	maxw = max(maxw,w++);
      }

  initmx((*M),A->m+1,maxw+1);

  for(i=1;i<= A->m;i++) for(j=0; j< A->n; j++) {
      J = A->IA[i-1][j];
      if (J != 0) mx((*M),i,J) = A->A[i-1][j];
    }
}

