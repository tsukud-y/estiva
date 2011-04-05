#include <stdio.h>
#include "estiva/vec.h"
#include "estiva/ary.h"
#include "estiva/std.h"

int estiva_fprintvec(double *b, char *name)
{
  FILE *fp;
  long i;
  fp = fopen(name,"w");
  if ( fp == NULL ) return 1;
  forall(0,i,dim1(b)) fprintf(fp,"%e\n",b[i]);
  fclose(fp);
  return 0;
}
