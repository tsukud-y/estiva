#include <stdio.h>
#include <stdlib.h>
#include "estiva/std.h"
#include "estiva/ary.h"

void estiva_swap(void *A, void *B, int size1, int size2)
{ static char *p3, *p1, *p2;static int i;
  p1=(char *)A;
  p2=(char *)B;

  ary1(p3,size1);
  
  if(size1 == size2 ) {
    forall(0,i,size1-1){ p3[i]=p2[i]; p2[i]=p1[i]; p1[i]=p3[i];}
  }
  else{ fprintf(stderr,"estiva_swap: size mismatch\n"); exit(1);}
}

void estiva_cp(void *A, void *B, int size1, int size2)
{ static int i; char *p1, *p2;
  p1=(char *)A;
  p2=(char *)B;

  if(size1 == size2) forall(0,i,size1-1) p2[i]=p1[i]; 
  else{ fprintf(stderr,"estiva_cp: size mismatch\n"); exit(1);}
}

