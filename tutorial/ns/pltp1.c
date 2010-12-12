#include "fem.h"

#include <stdio.h>
#include "estiva/mesh.h"
#include "estiva/foreach.h"
#include "estiva/ary.h"


static void pltmshx(FILE *fp, xyc *Z, nde *N, double *x)
{

  long e, a, b, c, dim1N, m;
  m = dimp2(N);
  dim1N = dim1(N);
  for(e=1;e<=dim1N;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y, x[2*m+a]);
    fprintf(fp,"%f %f %f\n",Z[b].x,Z[b].y, x[2*m+b]);
    fprintf(fp,"%f %f %f\n",Z[c].x,Z[c].y, x[2*m+c]);
    fprintf(fp,"%f %f %f\n",Z[a].x,Z[a].y, x[2*m+a]);
    fprintf(fp,"\n\n");
  }
}

void estiva_pltp1(double *x, xyc * Z, nde *N)
{
  FILE *pp;
  static double *u, *v;
  long e, m, dim1N, i;
  int n;
  double x0, y0;

  pp = popen("gnuplot","w");

  u = x;
  v = x + dimp2(N);

  for (i=1; i<=dimp2(N); i++) u[i]/=10.0;
  for (i=1; i<=dimp2(N); i++) v[i]/=10.0;

  dim1N = dim1(N);
  for (e=1; e<=dim1N; e++)
    foreach(n) &N[e].a, &N[e].b, &N[e].c, end {
      fprintf(pp,"set arrow from %f,%f to %f,%f\n",Z[n].x,Z[n].y,Z[n].x+u[n],Z[n].y+v[n]);
    }

  for (e=1; e<=dim1N; e++) {
    m = N[e].A;
    x0 = (Z[N[e].b].x + Z[N[e].c].x)/2.0, y0 = (Z[N[e].b].y + Z[N[e].c].y)/2.0;
    fprintf(pp,"set arrow from %f,%f to %f,%f\n",x0,y0,x0+u[m],y0+v[m]);

    m = N[e].B;
    x0 = (Z[N[e].c].x + Z[N[e].a].x)/2.0, y0 = (Z[N[e].c].y + Z[N[e].a].y)/2.0;
    fprintf(pp,"set arrow from %f,%f to %f,%f\n",x0,y0,x0+u[m],y0+v[m]);

    m = N[e].C;
    x0 = (Z[N[e].a].x + Z[N[e].b].x)/2.0, y0 = (Z[N[e].a].y + Z[N[e].b].y)/2.0;
    fprintf(pp,"set arrow from %f,%f to %f,%f\n",x0,y0,x0+u[m],y0+v[m]);
  }

  fprintf(pp,"splot '-' title \"\" with lines\n");
  pltmshx(pp,Z,N,x);
  fprintf(pp,"e\n");
  fflush(pp);
}
