#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/std.h"

#define y(x) (*(long *)f(x))

void estiva_forgammap1(long *x)
{
  if ( &y(x) == NULL ) Rnew(x, long); 
  if ( &y(x) == NULL ) abort();
  y(x) = 1;
}

int estiva_forgammap1_loop(long *x, char *NAME, xyc *Z)
{
  long  n;
  n = dim1(Z);

  while ( y(x) <= n ) {
    if (Z[y(x)].label && !strcmp(Z[y(x)].label,NAME) ) {
      *x = y(x);
      y(x)++;
      return 1;
    }
    else {
      y(x)++;
    }
  }
  
  Rdestroy(x);
  return 0;
}
