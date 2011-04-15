#include <stdio.h>
#include <stdlib.h>
#include "estiva/mesh.h"
#include "estiva/que.h"
#include "estiva/op.h"


static int eq(double a, double b)
{
  double eps = 0.0009, c;
  c = a-b;
  if ( c < 0.0 ) c = -c;
  if ( c <= eps ) return 1;
  return 0;
}



void estiva_rectmesh(xyc **Zp, nde **Np)
{
  static xyc *Z;
  static nde *N;

  static que *q;
  double x, y, h;
  long i, j;
  
  initq(q);

  h = 0.125;
  if ( defop("-n") ) h = 1.0/atof(getop("-n"));
  if ( defop("-h") ) h = atof(getop("-h"));

  for ( i = 0, x = 0.0; i < 1.0/h+h/2.0; i++, x += h)
    for ( j = 0, y = 0.0; j < 1.0/h+h/2.0; j++, y += h) {
      if ( x == 0.0 || eq(x,1.0) || y == 0.0 || eq(y,1.0) ) {
	if ( eq(y,1.0) && x != 0.0 && !eq(x,1.0) )
	  pushxyc(q,x,y,"north");
	else if ( y == 0.0 && x != 0.0 && !eq(x,1.0) )
	  pushxyc(q,x,y,"south");
	else if ( x == 0.0 && y != 0.0 && !eq(y,1.0) )
	  pushxyc(q,x,y,"west");
	else if ( eq(x,1.0) && y != 0.0 && !eq(y,1.0) )
	  pushxyc(q,x,y,"east");
	else 
	  pushxyc(q,x,y,"zero");
      }
      else 
	pushxyc(q,x,y,NULL);
    }
  pushxyc(q,h/2.0,h/2.0,NULL);
  pushxyc(q,1.0-h/2.0,1.0-h/2.0,NULL);
  
  genmesh(q,Z,N);
  p2(Z,N);
  
  *Zp = Z; *Np = N;
}

