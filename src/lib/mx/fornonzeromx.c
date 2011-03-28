#include "estiva/std.h"
#include "estiva/mx.h"

#define ii static_bind(long,Ip)
#define jj static_bind(long,Jp)

void estiva_fornonzeromx(MX *A){
  mx(A,1,1) = mx(A,1,1);
}

int estiva_fornonzeromx_loop(MX *A, long *Ip, long *Jp){
  long i, j;

  for (i = ii; i < A->m; i++, ii++, jj=0) 
    for (j = jj; j< A->n; j++) 
      if (A->IA[i][j] != 0.0) {
	*Ip = i+1;
	*Jp = A->IA[i][j];
	jj = j+1;
	return 1;
      }


  static_free(Ip);
  static_free(Jp);
  return 0;
}  
