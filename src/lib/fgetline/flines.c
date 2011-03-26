#include <stdio.h>
#include "estiva/fgetline.h"


long estiva_flines(void *vfp)
{
  FILE *fp;
  long i;
  fp = vfp;

  for ( i=0; !feof(fp); i++ ) fgetline(fp);
  rewind(fp);
  return i - 1;
}

