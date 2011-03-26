#include <stdio.h>
#include "estiva/fgetline.h"

long estiva_fsize(void *vfp)
{
  FILE *fp;
  long i=0;
  fp = vfp;
  while( EOF != fgetc(fp)) i++;
  return i;
}
