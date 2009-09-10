#include <stdio.h>
#include <stdlib.h>
#include "estiva.h"
#include "confary.h"
#include "matprop.h"
#include "solver.h"

static int *LU(void *A)
{ if(dim2(A) == dim1(A))
    if(siz(A) == sizeof(float)){ float nrm,**a,aik,*apk,*ai,apkk;
#include "solver/LU.c"
    }else if(siz(A) == sizeof(double)){ double nrm,**a,aik,*apk,*ai,apkk;
#include "solver/LU.c"
    }
  return NULL;
}

int solver_gauss(void *A, void *B)
{ void *X;

  if(dim1(A) != dim1(B)) return 0;
  { static int **Y;
    if(-1 == dim2(B)){
      ary2(Y,2,1);
      if(NULL == Y) return 0;
      cp(B,Y[1]); cp(Y,X);
    }else cp(B,X);
  }
  if(siz(A) == sizeof(float) && siz(B) == sizeof(float)){
    float **a,*api,s,**x,*b,bpi;
#include "solver/fb.c"
  }else if(siz(A) == sizeof(float) && siz(B) == sizeof(double)){
    float **a,*api;double s,**x,*b,bpi;
#include "solver/fb.c"
  }else if(siz(A) == sizeof(double) && siz(B) == sizeof(float)){
    double **a,*api,s;float **x,*b,bpi;
#include "solver/fb.c"
  }else	if(siz(A) == sizeof(double) && siz(B) == sizeof(double)){
    double **a,*api,s,**x,*b,bpi;
#include "solver/fb.c"
  }else{ fprintf(stderr,"solver_gauss: not supported\n"); exit(1);}
  return 0;
}

int solver_inv(void *A)
{ if(dim2(A) == dim1(A)) 
    if(siz(A) == sizeof(float)){ static float **a,**e;
#include "solver/inv.c"    
    }else if(siz(A) == sizeof(double)){ static double **a,**e;
#include "solver/inv.c"    
    }else{ fprintf(stderr,"solver_inv: not supported\n"); exit(1);}
  return 0;
}
