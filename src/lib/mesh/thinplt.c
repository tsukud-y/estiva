#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/mesh.h"
#include "estiva/ary.h"
#include "estiva/op.h"


static double fplane(double x, double y, double z,
                     double x0,double y0,double z0,
                     double x1,double y1,double z1,
                     double x2,double y2,double z2)
{ double xn,yn,zn;

  xn = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
  yn = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
  zn = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);

  return xn*x+yn*y+zn*z-(xn*x0+yn*y0+zn*z0);
}


static int notLawson(double x, double y, xyc *Z, nde *N)
{
  int e,n;
  n = dim1(N);

  for(e=1;e<=n;e++){
    int a,b,c; double x0,y0,x1,y1,x2,y2;
    a = N[e].a; b = N[e].b; c = N[e].c; 
    x0 = Z[a].x; x1 = Z[b].x; x2 = Z[c].x;
    y0 = Z[a].y; y1 = Z[b].y; y2 = Z[c].y;

    if     (fplane(x,y,0.0, x0,y0,1.0, x1,y1,0.0, x2,y2,0.0)>0.0) continue;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,1.0, x2,y2,0.0)>0.0) continue;
    else if(fplane(x,y,0.0, x0,y0,0.0, x1,y1,0.0, x2,y2,1.0)>0.0) continue;
    else return e;
  }
  return  0;
}


static double alpha(int e, double x, double y, xyc *Z, nde *N)
{
  double xa, xb, xc, ya, yb, yc, delta;
  
  xa = Z[N[e].a].x;
  xb = Z[N[e].b].x;
  xc = Z[N[e].c].x;
  ya = Z[N[e].a].y;
  yb = Z[N[e].b].y;
  yc = Z[N[e].c].y;

  delta = (xa-xc)*(yb-yc)-(ya-yc)*(xb-xc);

  if ( delta == 0.0) return 1.0;

  return ((x-xc)*(yb-yc)-(y-yc)*(xb-xc))/delta;
}

static double beta(int e, double x, double y, xyc *Z, nde *N)
{
  double xa, xb, xc, ya, yb, yc, delta;
  
  xa = Z[N[e].a].x;
  xb = Z[N[e].b].x;
  xc = Z[N[e].c].x;
  ya = Z[N[e].a].y;
  yb = Z[N[e].b].y;
  yc = Z[N[e].c].y;

  delta = (xb-xc)*(ya-yc)-(yb-yc)*(xa-xc);

  if ( delta == 0.0 ) return 1.0;

  return ((x-xc)*(ya-yc)-(y-yc)*(xa-xc))/delta;
}

static char *num2str(long n){
  long i;
  static char ret[100];
  sprintf(ret,"00000000%ld",n);
  i = strlen( ret );
  
  return &ret[i-8];
}


void estiva_thinplt(double *u, xyc *Z, nde *N)
{
  FILE *pp;
  static char filename[100];
  static long arrow=1, sequence=0;
  
  int find, i, j; double x, y, a, b, *v, scale = 0.001;
  v = u + dimp2(N);

  if ( defop("-plotscale") ) scale = atof(getop("-plotscale"));
  
  system("mkdir -p anime");
  sprintf(filename,"anime/%s.gnuplot",num2str(sequence++));
  pp = fopen(filename,"w");

  for (x=0.0,j=0; j<32; x+=0.03125,j++) 
    for (y=0.0, i=0; i<32; y+=0.03125, i++){
      find = notLawson(x,y,Z,N);
      if ( find == 0 ) continue;      
      a = alpha(find,x,y,Z,N);
      b = beta(find, x,y,Z,N);
      fprintf(pp,"set arrow %ld from %f,%f to %f,%f\n",
	      arrow++,x,y,
	      x+scale*(a*u[N[find].a]+b*u[N[find].b]+(1.0-a-b)*u[N[find].c]),
	      y+scale*(a*v[N[find].a]+b*v[N[find].b]+(1.0-a-b)*v[N[find].c]));
    }
  arrow = 1;
  fprintf(pp,"plot '-' title \"\" with lines\n");
  fprintf(pp,"0.0 0.0\n");
  fprintf(pp,"1.0 0.0\n");
  fprintf(pp,"1.0 1.0\n");
  fprintf(pp,"0.0 1.0\n");
  fprintf(pp,"\n");
  fprintf(pp,"e\n");
  fprintf(pp,"\n");
  fflush(pp);
  fclose(pp);
}

