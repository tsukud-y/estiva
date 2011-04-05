#include <stdio.h>
#include "estiva/mx.h"

int estiva_fprintmx(void *Apointer, char *name)
{
  FILE *fp;
  MX *A;
  long i, j;
  A = Apointer;
  fp = fopen(name,"w");
  if ( fp == NULL ) return 1;
  fornonzeromx(A,i,j) fprintf(fp,"%ld %ld %e\n",i,j,mx(A,i,j));
  fclose(fp);
  return 0;
}
