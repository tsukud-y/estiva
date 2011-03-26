#include <stdio.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"

static long n, i;
static char *NAME;
static xyc *Z;

void estiva_forgammap1(char *NAMEp, xyc *Zp)
{
  Z = Zp;
  n = dim1(Z);
  i = 1;
  NAME = NAMEp;
}


int estiva_forgammap1_loop(long *ip)
{
  while ( i <= n ) {
    if (Z[i].label && !strcmp(Z[i].label,NAME) ) {
      *ip = i;
      i++;
      return 1;
    }
    else {
      i++;
    }
  }
  return 0;
}
