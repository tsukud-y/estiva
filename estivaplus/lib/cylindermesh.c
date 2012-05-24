#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/mesh.h"
#include "estiva/que.h"
#include "estiva/op.h"

void estiva_setmesh(xyc *Z, nde *N);

#define fortheta(theta,T)					\
  for ( theta = 0.0; theta <  2.0*M_PI-(1.0*M_PI)/T; theta += (2.0*M_PI)/T)

void cylindermesh(xyc **Zp, nde **Np)
{
  static xyc *Z;
  static nde *N;

  static que *q;
  double x, y;
  double theta, r;
  double T = 32.0;

  initq(q);
  
  fortheta(theta,T) pushxyc(q, cos(theta), sin(theta),"cylinder");

  for ( r = 1.1; r <= 1.5; r += 0.1)
    fortheta(theta,T) pushxyc(q, r*cos(theta), r*sin(theta),NULL);

  for ( r = 2.0; r <= 7.0; r += 0.5)
    fortheta(theta,T) pushxyc(q, r*cos(theta), r*sin(theta),NULL);
  
  for ( r = 8.0; r <= 10.0; r += 1.0)
    fortheta(theta,T)  {
      x = r*cos(theta); y = r*sin(theta);
      if (-7.41 < x && x < 7.41 && -7.41 < y && y < 7.41) 
	pushxyc(q, x, y, NULL);
    }  
  
  pushxyc(q, -7.5,  7.5, "zero");
  pushxyc(q, -7.5, -7.5, "zero");
  
  for ( x = -6.5; x <= 7.5; x += 1.0 ) {
    pushxyc(q, x,  7.5, "north");
    pushxyc(q, x, -7.5, "south");
  }
  for ( y = -6.5; y <= 6.5; y += 1.0 ) {
    pushxyc(q, -7.5, y, "west");
    pushxyc(q,  7.5, y, "east");
  }

  genmesh(q,Z,N);

  p2(Z,N);
  estiva_setmesh(Z,N);

  *Zp = Z; *Np = N;
}
