#include "estiva/mx.h"


void estiva_fornonzeromx(MX *A){
  mx(A,1,1) = mx(A,1,1);
}


int estiva_fornonzeromx_loop(MX *A, long *Ip, long *Jp){
  static long ii=0, jj=0;
  long i, j;

  for (i = ii; i < A->m; i++, ii++, jj=0) 
    for (j = jj; j< A->n; j++) 
      if (A->IA[i][j] != 0.0) {
	*Ip = i+1;
	*Jp = A->IA[i][j];
	jj = j+1;
	return 1;
      }

  ii=0, jj=0;
  return 0;
}  
