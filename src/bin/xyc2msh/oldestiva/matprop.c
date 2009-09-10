#include <stdio.h>
#include "estiva.h"
#include "confary.h"
#include "matprop.h"

int matprop_halfbw(void *A)
{ 
  if(dim2(A) != dim1(A)||dim2(A) < 1||dim1(A) < 1||siz(A) < 1) return 0;
  
  if(siz(A) == sizeof(float)){
    static float **a; 
    a = (float **)A;
#include "matprop/halfbw.c"
  }
  if(siz(A) == sizeof(double)){
    static double **a; 
    a = (double **)A;
#include "matprop/halfbw.c"
  }
  return 0;
}
