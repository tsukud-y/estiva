#include <stdio.h>
#include "estiva/fgetline.h"


long estiva_flines(FILE *fp)
{
  long i;
  for ( i=0; !feof(fp); i++ ) fgetline(fp);
  rewind(fp);
  return i - 1;
}

