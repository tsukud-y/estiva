#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/que.h"


void estiva_pushxyc(void *p, double x, double y, char *label)
{
  xyc e;
  que **qp;
  qp = p;
  e.x = x; e.y = y; e.label = label;
  push(*qp,e);
}
