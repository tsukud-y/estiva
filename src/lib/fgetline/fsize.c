#include <stdio.h>

long estiva_fsize(FILE *fp)
{
  long i=0;
  while( EOF != fgetc(fp)) i++;
  return i;
}
