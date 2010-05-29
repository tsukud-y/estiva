#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/que.h"


char *estiva_fgetline(FILE *fp)
{
  static char *buf, *ret;
  static que *q;
  unsigned long i, j, n=80;
  int c; 

  ary1(buf,n+1);
  initq(q); while ( NULL != pop(q,c) );
  
  for ( i=0; EOF != (c = fgetc(fp)) && c != '\n'; )
    if ( i < n )  buf[i++] = c;
    else          enq(q,c);
  buf[i] = '\n';
  buf[i+1] = 0;
  j = 0; forq(q,c) j++;
  
  ary1(ret,n+j+2);
  ret[0] = 0;
  strcat(ret,buf);
  
  forq(q,c) ret[i++] = c;
  ret[i] = '\n';
  ret[i+1] = 0;

  return ret;
}
