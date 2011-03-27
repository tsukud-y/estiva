#include <stdio.h>
#include "estiva/op.h"


void *estiva_ofp(void)
{
  FILE *fp;

  if (NULL != (fp=fopen(getop("-o"),"w")) )
    return fp;
  return stdout;
}
