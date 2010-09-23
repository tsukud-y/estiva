#include "fem.h"

#include <stdlib.h>
#include "estiva/que.h"
#include "estiva/op.h"


void estiva_rectmesh(xyc **Zp, nde **Np)
{
  static xyc *Z;
  static nde *N;

  static que *q;
  double x, y, h;
  
  initq(q);

  h = 0.125;
  if ( defop("-h") ) h = atof(getop("-h"));
  for ( x = 0.0; x <= 1.0; x += h)
    for ( y = 0.0; y <= 1.0; y += h) {
      
      if ( x == 0.0 || x == 1.0 || y == 0.0 || y == 1.0 ) {
	if ( y == 1.0 && x != 0.0 && x != 1.0 )
	  pushxyc(q,x,y,"gamma");
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

