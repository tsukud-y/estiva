#include "estivaplus.h"

static xyc *Z_easymesh;
static nde *N_easymesh;
static que *q;

void estiva_putnode(double x, double y, char *label)
{
  static long n=0;
  if ( n == 0 ) {
    initq(q);
  }
  n++;
  pushxyc(q,x,y,label);
}

void GenerateMesh(void)
{
  genmesh(q,Z_easymesh,N_easymesh);
  p2(Z_easymesh,N_easymesh);
  estiva_setmesh(Z_easymesh,N_easymesh);  
}

